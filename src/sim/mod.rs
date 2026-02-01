use macroquad::prelude::*;
use rand09::rngs::StdRng;
use rand09::{Rng, SeedableRng};

use crate::ui::Controls;

#[derive(Clone, Copy)]
pub struct Cell {
    pub temp: f32,
    pub oil: f32,
    pub vel: Vec2,
}

#[derive(Clone)]
pub struct Particle {
    pub p: Vec2,
    pub life: f32,
}

#[derive(Clone, Copy)]
pub struct BoomSeg {
    pub a: Vec2,
    pub b: Vec2,
    pub block: f32, // 0..1, how strongly it blocks crossing
}

#[derive(Clone, Copy)]
pub struct Skimmer {
    pub p: Vec2,
    pub r: f32,
    pub rate: f32, // higher = removes more per second
}

#[derive(Clone, Copy)]
pub struct Dispersant {
    pub p: Vec2,
    pub r: f32,
    pub strength: f32, // higher = more decay + more diffusion
}

pub struct World {
    pub w: usize,
    pub h: usize,
    pub cell: f32,

    pub grid: Vec<Cell>,
    pub particles: Vec<Particle>,

    pub booms: Vec<BoomSeg>,
    pub skimmers: Vec<Skimmer>,
    pub dispersants: Vec<Dispersant>,

    rng: StdRng,
    t: f32,
}

pub struct Metrics {
    pub oil_mass: f32,
    pub oil_max: f32,
}

impl World {
    pub fn new(w: usize, h: usize) -> Self {
        let mut rng = StdRng::seed_from_u64(12345);

        let cell = 7.0;
        let mut grid = vec![
            Cell { temp: 0.08, oil: 0.0, vel: vec2(0.0, 0.0) };
            w * h
        ];

        for y in 0..h {
            for x in 0..w {
                let fx = x as f32 / (w as f32 - 1.0);
                let fy = y as f32 / (h as f32 - 1.0);
                let base = 0.06 + 0.06 * (1.0 - fy) + 0.02 * (fx * (1.0 - fx));
                grid[y * w + x].temp = base;
            }
        }

        let mut particles = Vec::new();
        for _ in 0..5500 {
            particles.push(Particle {
                p: vec2(
                    rng.random_range(0.0..(w as f32 * cell)),
                    rng.random_range(0.0..(h as f32 * cell)),
                ),
                life: rng.random_range(0.3..1.0),
            });
        }

        Self {
            w,
            h,
            cell,
            grid,
            particles,
            booms: Vec::new(),
            skimmers: Vec::new(),
            dispersants: Vec::new(),
            rng,
            t: 0.0,
        }
    }

    pub fn vel_at(&self, p: Vec2) -> Vec2 {
        sample_vel_grid(&self.grid, self.w, self.h, self.cell, p)
    }

    pub fn metrics(&self) -> Metrics {
        let mut mass = 0.0;
        let mut mx: f32 = 0.0;
        for c in &self.grid {
            mass += c.oil;
            mx = mx.max(c.oil);
        }
        Metrics { oil_mass: mass, oil_max: mx }
    }

    pub fn step(&mut self, dt: f32, c: &Controls) {
        self.t += dt;

        self.update_velocity(dt, c);
        self.advect_scalars(dt, c);
        self.diffuse_temp(dt, c);
        self.diffuse_oil(dt, c);
        self.apply_response(dt, c);
        self.advect_particles(dt);
        self.weather_oil(dt, c);
        self.decay_temp(dt);
    }

    fn update_velocity(&mut self, dt: f32, c: &Controls) {
        let wind = vec2(
            c.wind_dir_deg.to_radians().cos(),
            c.wind_dir_deg.to_radians().sin(),
        ) * (120.0 * c.wind_strength);

        for y in 0..self.h {
            for x in 0..self.w {
                let i = y * self.w + x;

                let px = x as f32 / self.w as f32;
                let py = y as f32 / self.h as f32;

                let a = (px * 7.0 + self.t * 0.35).sin() * (py * 6.0 + self.t * 0.22).cos();
                let b = (py * 8.0 + self.t * 0.31).sin() * (px * 5.0 + self.t * 0.27).cos();
                let swirl = vec2(-b, a) * 80.0;

                let v = wind + swirl;
                self.grid[i].vel = self.grid[i].vel.lerp(v, (dt * 2.6).clamp(0.0, 1.0));
            }
        }
    }

    fn advect_scalars(&mut self, dt: f32, c: &Controls) {
        let a = c.advection.clamp(0.0, 2.0);
        if a <= 0.0 { return; }

        let prev = self.grid.clone();
        let mut next = self.grid.clone();

        let w = self.w;
        let h = self.h;
        let cell = self.cell;

        for y in 0..h {
            for x in 0..w {
                let i = y * w + x;

                let center = vec2((x as f32 + 0.5) * cell, (y as f32 + 0.5) * cell);
                let v = self.grid[i].vel;

                let back = center - v * (dt * a);

                let (t, oil) = sample_temp_oil(&prev, w, h, cell, back);

                next[i].temp = t;
                next[i].oil = oil;
                next[i].vel = self.grid[i].vel;
            }
        }

        // apply boom blocking after advection (cheap but works)
        self.grid = next;
        self.apply_booms(dt);
    }

    fn apply_booms(&mut self, dt: f32) {
        if self.booms.is_empty() { return; }

        // crude but effective:
        // for each boom, we damp oil transport across the boom by reducing oil in cells "crossing" direction.
        // we approximate by projecting cell centers and applying a local barrier.
        let wpx = self.w as f32 * self.cell;
        let hpx = self.h as f32 * self.cell;

        let mut next = self.grid.clone();

        for boom in &self.booms {
            let a = boom.a;
            let b = boom.b;
            let ab = b - a;
            let len = ab.length();
            if len < 1.0 { continue; }

            let n = vec2(-ab.y, ab.x) / len; // normal
            let minx = a.x.min(b.x) - 60.0;
            let maxx = a.x.max(b.x) + 60.0;
            let miny = a.y.min(b.y) - 60.0;
            let maxy = a.y.max(b.y) + 60.0;

            let minx = minx.clamp(0.0, wpx);
            let maxx = maxx.clamp(0.0, wpx);
            let miny = miny.clamp(0.0, hpx);
            let maxy = maxy.clamp(0.0, hpx);

            let x0 = (minx / self.cell) as i32;
            let x1 = (maxx / self.cell) as i32;
            let y0 = (miny / self.cell) as i32;
            let y1 = (maxy / self.cell) as i32;

            for yy in y0..=y1 {
                for xx in x0..=x1 {
                    if xx < 0 || yy < 0 || xx >= self.w as i32 || yy >= self.h as i32 { continue; }
                    let i = yy as usize * self.w + xx as usize;

                    let p = vec2((xx as f32 + 0.5) * self.cell, (yy as f32 + 0.5) * self.cell);

                    let d = signed_dist_to_segment(p, a, b, n);
                    if d.abs() > 18.0 { continue; }

                    // if local velocity tries to cross the boom, suppress oil near it
                    let v = self.grid[i].vel;
                    let cross = v.dot(n).abs();

                    let block = (boom.block * (cross / 140.0).clamp(0.0, 1.0)).clamp(0.0, 0.98);
                    let keep = 1.0 - block * (dt * 6.0).clamp(0.0, 1.0);

                    next[i].oil *= keep;
                }
            }
        }

        for i in 0..next.len() {
            next[i].vel = self.grid[i].vel;
            next[i].temp = self.grid[i].temp;
        }
        self.grid = next;
    }

    fn diffuse_temp(&mut self, dt: f32, c: &Controls) {
        let k = (c.diffusion * dt).clamp(0.0, 0.5);
        if k <= 0.0 { return; }
        let mut next = self.grid.clone();

        for y in 0..self.h {
            for x in 0..self.w {
                let i = y * self.w + x;

                let mut sum = self.grid[i].temp;
                let mut n = 1.0;

                for (dx, dy) in [(-1, 0), (1, 0), (0, -1), (0, 1)] {
                    let nx = x as i32 + dx;
                    let ny = y as i32 + dy;
                    if nx >= 0 && ny >= 0 && (nx as usize) < self.w && (ny as usize) < self.h {
                        let j = ny as usize * self.w + nx as usize;
                        sum += self.grid[j].temp;
                        n += 1.0;
                    }
                }

                let avg = sum / n;
                next[i].temp += (avg - self.grid[i].temp) * k;
            }
        }

        for i in 0..next.len() {
            next[i].vel = self.grid[i].vel;
            next[i].oil = self.grid[i].oil;
        }
        self.grid = next;
    }

    fn diffuse_oil(&mut self, dt: f32, c: &Controls) {
        let k = (c.oil_diffusion * dt).clamp(0.0, 0.5);
        if k <= 0.0 && self.dispersants.is_empty() { return; }

        let mut next = self.grid.clone();

        for y in 0..self.h {
            for x in 0..self.w {
                let i = y * self.w + x;

                // dispersant locally increases diffusion
                let p = vec2((x as f32 + 0.5) * self.cell, (y as f32 + 0.5) * self.cell);
                let mut dk = k;
                for d in &self.dispersants {
                    let r = d.r.max(1.0);
                    let t = 1.0 - (p.distance(d.p) / r).clamp(0.0, 1.0);
                    dk = (dk + dk * d.strength * 0.75 * t).clamp(0.0, 0.5);
                }
                if dk <= 0.0 { continue; }

                let mut sum = self.grid[i].oil;
                let mut n = 1.0;

                for (dx, dy) in [(-1, 0), (1, 0), (0, -1), (0, 1)] {
                    let nx = x as i32 + dx;
                    let ny = y as i32 + dy;
                    if nx >= 0 && ny >= 0 && (nx as usize) < self.w && (ny as usize) < self.h {
                        let j = ny as usize * self.w + nx as usize;
                        sum += self.grid[j].oil;
                        n += 1.0;
                    }
                }

                let avg = sum / n;
                next[i].oil += (avg - self.grid[i].oil) * dk;
                next[i].oil = next[i].oil.max(0.0);
            }
        }

        for i in 0..next.len() {
            next[i].vel = self.grid[i].vel;
            next[i].temp = self.grid[i].temp;
        }
        self.grid = next;
    }

    fn apply_response(&mut self, dt: f32, c: &Controls) {
        // skimmers remove oil in radius
        if !self.skimmers.is_empty() {
            for sk in &self.skimmers {
                let r = sk.r.max(1.0);
                let r2 = r * r;
                let rc = (r / self.cell).ceil() as i32 + 2;

                let cx = (sk.p.x / self.cell) as i32;
                let cy = (sk.p.y / self.cell) as i32;

                for dy in -rc..=rc {
                    for dx in -rc..=rc {
                        let x = cx + dx;
                        let y = cy + dy;
                        if x < 0 || y < 0 || x >= self.w as i32 || y >= self.h as i32 { continue; }
                        let i = y as usize * self.w + x as usize;

                        let p = vec2((x as f32 + 0.5) * self.cell, (y as f32 + 0.5) * self.cell);
                        let d2 = p.distance_squared(sk.p);
                        if d2 > r2 { continue; }

                        let fall = 1.0 - (d2 / r2).sqrt();
                        let rm = (sk.rate * fall * dt).clamp(0.0, 0.8);
                        self.grid[i].oil *= 1.0 - rm;
                    }
                }
            }
        }

        // optional: continuous source at mouse? handled by UI splat; but we also support "source rate" via a hotkey later if you want
        let _ = c.oil_source_rate;
    }

    fn weather_oil(&mut self, dt: f32, c: &Controls) {
        // base weathering
        for y in 0..self.h {
            for x in 0..self.w {
                let i = y * self.w + x;

                let p = vec2((x as f32 + 0.5) * self.cell, (y as f32 + 0.5) * self.cell);

                // dispersant locally increases decay
                let mut decay = c.oil_decay;
                for d in &self.dispersants {
                    let r = d.r.max(1.0);
                    let t = 1.0 - (p.distance(d.p) / r).clamp(0.0, 1.0);
                    decay += c.oil_decay * d.strength * 1.25 * t;
                }

                let keep = (1.0 - decay * dt).clamp(0.0, 1.0);
                self.grid[i].oil *= keep;
            }
        }
    }

    fn decay_temp(&mut self, dt: f32) {
        let td = (1.0 - dt * 0.02).clamp(0.94, 1.0);
        for c in &mut self.grid {
            c.temp *= td;
        }
    }

    fn advect_particles(&mut self, dt: f32) {
        let maxx = self.w as f32 * self.cell;
        let maxy = self.h as f32 * self.cell;

        let w = self.w;
        let h = self.h;
        let cell = self.cell;

        let grid = &self.grid;
        let rng = &mut self.rng;
        let particles = &mut self.particles;

        for pt in particles.iter_mut() {
            let v = sample_vel_grid(grid, w, h, cell, pt.p);
            pt.p += v * dt;

            if pt.p.x < 0.0 { pt.p.x += maxx; }
            if pt.p.x >= maxx { pt.p.x -= maxx; }
            if pt.p.y < 0.0 { pt.p.y += maxy; }
            if pt.p.y >= maxy { pt.p.y -= maxy; }

            pt.life += dt * 0.25;
            if pt.life > 1.0 {
                pt.life = rng.random_range(0.2..0.6);
            }
        }
    }

    pub fn splat_temp(&mut self, pos: Vec2, r: f32, amt: f32) {
        self.splat(pos, r, amt, 0.0);
    }

    pub fn splat_oil(&mut self, pos: Vec2, r: f32, amt: f32) {
        self.splat(pos, r, 0.0, amt);
    }

    fn splat(&mut self, pos: Vec2, r: f32, add_t: f32, add_oil: f32) {
        let r2 = r * r;

        let cx = (pos.x / self.cell) as i32;
        let cy = (pos.y / self.cell) as i32;
        let rc = (r / self.cell).ceil() as i32 + 1;

        for dy in -rc..=rc {
            for dx in -rc..=rc {
                let x = cx + dx;
                let y = cy + dy;

                if x < 0 || y < 0 || x >= self.w as i32 || y >= self.h as i32 {
                    continue;
                }

                let wx = (x as f32 + 0.5) * self.cell;
                let wy = (y as f32 + 0.5) * self.cell;
                let d2 = (wx - pos.x).powi(2) + (wy - pos.y).powi(2);
                if d2 > r2 { continue; }

                let f = 1.0 - (d2 / r2).sqrt();
                let i = y as usize * self.w + x as usize;

                self.grid[i].temp = (self.grid[i].temp + add_t * f).clamp(0.0, 1.0);
                self.grid[i].oil = (self.grid[i].oil + add_oil * f).max(0.0);
            }
        }
    }

    pub fn erase_response_near(&mut self, p: Vec2, r: f32) {
        let r2 = r * r;
        self.skimmers.retain(|s| s.p.distance_squared(p) > r2);
        self.dispersants.retain(|d| d.p.distance_squared(p) > r2);

        // for booms: distance to segment
        self.booms.retain(|b| dist_point_to_segment2(p, b.a, b.b) > r2);
    }
}

fn sample_vel_grid(grid: &[Cell], w: usize, h: usize, cell: f32, p: Vec2) -> Vec2 {
    let fx = (p.x / cell).clamp(0.0, w as f32 - 1.001);
    let fy = (p.y / cell).clamp(0.0, h as f32 - 1.001);

    let x0 = fx.floor() as usize;
    let y0 = fy.floor() as usize;
    let x1 = (x0 + 1).min(w - 1);
    let y1 = (y0 + 1).min(h - 1);

    let tx = fx - x0 as f32;
    let ty = fy - y0 as f32;

    let v00 = grid[y0 * w + x0].vel;
    let v10 = grid[y0 * w + x1].vel;
    let v01 = grid[y1 * w + x0].vel;
    let v11 = grid[y1 * w + x1].vel;

    v00.lerp(v10, tx).lerp(v01.lerp(v11, tx), ty)
}

fn sample_temp_oil(grid: &[Cell], w: usize, h: usize, cell: f32, p: Vec2) -> (f32, f32) {
    let fx = (p.x / cell).clamp(0.0, w as f32 - 1.001);
    let fy = (p.y / cell).clamp(0.0, h as f32 - 1.001);

    let x0 = fx.floor() as usize;
    let y0 = fy.floor() as usize;
    let x1 = (x0 + 1).min(w - 1);
    let y1 = (y0 + 1).min(h - 1);

    let tx = fx - x0 as f32;
    let ty = fy - y0 as f32;

    let c00 = grid[y0 * w + x0];
    let c10 = grid[y0 * w + x1];
    let c01 = grid[y1 * w + x0];
    let c11 = grid[y1 * w + x1];

    let t0 = c00.temp + (c10.temp - c00.temp) * tx;
    let t1 = c01.temp + (c11.temp - c01.temp) * tx;
    let temp = t0 + (t1 - t0) * ty;

    let o0 = c00.oil + (c10.oil - c00.oil) * tx;
    let o1 = c01.oil + (c11.oil - c01.oil) * tx;
    let oil = o0 + (o1 - o0) * ty;

    (temp, oil)
}

fn signed_dist_to_segment(p: Vec2, a: Vec2, _b: Vec2, n: Vec2) -> f32 {
    // signed distance using segment normal, only meaningful near the segment
    let ap = p - a;
    ap.dot(n)
}

fn dist_point_to_segment2(p: Vec2, a: Vec2, b: Vec2) -> f32 {
    let ab = b - a;
    let t = ((p - a).dot(ab) / ab.dot(ab)).clamp(0.0, 1.0);
    let q = a + ab * t;
    p.distance_squared(q)
}

use macroquad::prelude::*;
use crate::sim::World;
use crate::ui::Controls;

pub fn draw_world(world: &World, c: &Controls, cam: &Camera2D) {
    set_camera(cam);

    if c.show_temp || c.show_oil {
        draw_heatmaps(world, c);
    }

    if c.show_streamlines {
        draw_streamlines(world, c);
    }

    if c.show_vectors {
        draw_vectors(world);
    }

    if c.show_particles {
        draw_particles(world);
    }

    if c.show_response {
        draw_response(world);
    }

    let w = world.w as f32 * world.cell;
    let h = world.h as f32 * world.cell;
    draw_rectangle_lines(0.0, 0.0, w, h, 2.0, GRAY);

    set_default_camera();
}

fn draw_heatmaps(world: &World, c: &Controls) {
    for y in 0..world.h {
        for x in 0..world.w {
            let i = y * world.w + x;
            let cell = world.grid[i];

            let mut col = Color::new(0.02, 0.02, 0.03, 1.0);

            if c.show_temp {
                col = blend(col, temp_col(cell.temp), 0.92);
            }

            if c.show_oil {
                let oil = (cell.oil / c.oil_max.max(0.001)).clamp(0.0, 1.0);
                col = blend(col, oil_col(oil), 0.95);
            }

            let px = x as f32 * world.cell;
            let py = y as f32 * world.cell;
            draw_rectangle(px, py, world.cell + 0.6, world.cell + 0.6, col);
        }
    }
}

fn oil_col(o: f32) -> Color {
    // dark slick w/ subtle iridescent edge
    let a = (o * 0.92).clamp(0.0, 0.92);
    let base = Color::new(0.02, 0.02, 0.03, a);
    let edge = Color::new(0.35, 0.55, 0.9, (o * 0.25).clamp(0.0, 0.25));
    blend(base, edge, (o * o).clamp(0.0, 1.0))
}

fn draw_response(world: &World) {
    for b in &world.booms {
        draw_line(b.a.x, b.a.y, b.b.x, b.b.y, 4.0, Color::new(1.0, 0.85, 0.2, 0.95));
        draw_circle(b.a.x, b.a.y, 3.0, Color::new(1.0, 0.85, 0.2, 0.95));
        draw_circle(b.b.x, b.b.y, 3.0, Color::new(1.0, 0.85, 0.2, 0.95));
    }

    for s in &world.skimmers {
        draw_circle_lines(s.p.x, s.p.y, s.r * 0.35, 2.0, Color::new(0.85, 0.95, 1.0, 0.85));
        draw_circle(s.p.x, s.p.y, 4.0, Color::new(0.85, 0.95, 1.0, 0.85));
    }

    for d in &world.dispersants {
        draw_circle_lines(d.p.x, d.p.y, d.r * 0.45, 2.0, Color::new(0.6, 1.0, 0.6, 0.55));
        draw_circle(d.p.x, d.p.y, 4.0, Color::new(0.6, 1.0, 0.6, 0.55));
    }
}

fn draw_streamlines(world: &World, c: &Controls) {
    let w = world.w as f32 * world.cell;
    let h = world.h as f32 * world.cell;

    let spacing = c.sl_spacing.max(10.0);
    let steps = c.sl_steps.max(3) as usize;
    let dt = c.sl_dt.clamp(0.001, 0.08);

    let mut y = spacing * 0.5;
    while y < h {
        let mut x = spacing * 0.5;
        while x < w {
            trace_streamline(world, vec2(x, y), steps, dt, w, h);
            x += spacing;
        }
        y += spacing;
    }
}

fn trace_streamline(world: &World, start: Vec2, steps: usize, dt: f32, w: f32, h: f32) {
    let mut p = start;
    let mut prev = p;

    for _ in 0..steps {
        let v1 = world.vel_at(p);
        let s1 = v1.length();
        if s1 < 6.0 { break; }

        let mid = p + v1 * (dt * 0.5);
        let v2 = world.vel_at(mid);
        let np = p + v2 * dt;

        if np.x < 0.0 || np.y < 0.0 || np.x >= w || np.y >= h { break; }

        let a = (s1 / 220.0).clamp(0.06, 0.42);
        draw_line(prev.x, prev.y, np.x, np.y, 1.2, Color::new(0.93, 0.93, 0.95, a));

        prev = np;
        p = np;
    }
}

fn draw_vectors(world: &World) {
    let step = 8usize;
    for y in (0..world.h).step_by(step) {
        for x in (0..world.w).step_by(step) {
            let i = y * world.w + x;
            let v = world.grid[i].vel;

            let p = vec2(
                (x as f32 + 0.5) * world.cell,
                (y as f32 + 0.5) * world.cell,
            );

            let mag = v.length();
            if mag < 1.0 { continue; }

            let dir = v / mag;
            let len = (mag / 140.0).clamp(0.25, 1.3) * (step as f32 * world.cell * 0.55);

            let q = p + dir * len;
            draw_line(p.x, p.y, q.x, q.y, 1.6, Color::new(0.92, 0.92, 0.92, 0.55));

            let n = vec2(-dir.y, dir.x);
            let a = q - dir * 6.0 + n * 3.0;
            let b = q - dir * 6.0 - n * 3.0;
            draw_triangle(q, a, b, Color::new(0.92, 0.92, 0.92, 0.55));
        }
    }
}

fn draw_particles(world: &World) {
    for pt in &world.particles {
        let a = (0.15 + 0.55 * pt.life).clamp(0.05, 0.8);
        draw_circle(pt.p.x, pt.p.y, 1.2, Color::new(0.9, 0.9, 0.95, a));
    }
}

fn temp_col(t: f32) -> Color {
    let t = t.clamp(0.0, 1.0);
    if t < 0.5 {
        let u = t / 0.5;
        lerp_color(Color::new(0.02, 0.12, 0.28, 1.0), Color::new(0.0, 0.75, 0.88, 1.0), u)
    } else {
        let u = (t - 0.5) / 0.5;
        lerp_color(Color::new(0.0, 0.75, 0.88, 1.0), Color::new(0.98, 0.55, 0.12, 1.0), u)
    }
}

fn lerp_color(a: Color, b: Color, t: f32) -> Color {
    Color::new(
        a.r + (b.r - a.r) * t,
        a.g + (b.g - a.g) * t,
        a.b + (b.b - a.b) * t,
        a.a + (b.a - a.a) * t,
    )
}

fn blend(base: Color, over: Color, k: f32) -> Color {
    let k = k.clamp(0.0, 1.0);
    Color::new(
        base.r + (over.r - base.r) * k * over.a,
        base.g + (over.g - base.g) * k * over.a,
        base.b + (over.b - base.b) * k * over.a,
        1.0,
    )
}

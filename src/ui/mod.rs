use crate::sim::{World, BoomSeg, Skimmer, Dispersant};
use egui_macroquad::egui;
use macroquad::prelude::*;

#[derive(Clone)]
pub struct Controls {
    pub paused: bool,
    pub sim_speed: f32,

    pub show_vectors: bool,
    pub show_particles: bool,
    pub show_temp: bool,
    pub show_oil: bool,
    pub show_streamlines: bool,
    pub show_response: bool,

    pub wind_dir_deg: f32,
    pub wind_strength: f32,

    pub advection: f32,
    pub diffusion: f32,

    pub oil_diffusion: f32,
    pub oil_decay: f32,
    pub oil_source_rate: f32,
    pub oil_max: f32,

    pub tool: ToolMode,
    pub brush_radius: f32,

    pub sl_spacing: f32,
    pub sl_steps: i32,
    pub sl_dt: f32,

    pub boom_blocks: f32,
    pub skimmer_rate: f32,
    pub dispersant_strength: f32,
}

#[derive(Clone, Copy, PartialEq)]
pub enum ToolMode {
    PaintHeat,
    SpillOil,
    PlaceBoom,
    PlaceSkimmer,
    PlaceDispersant,
    EraseResponse,
}

pub struct App {
    pub world: World,
    pub controls: Controls,
    pub camera: Camera2D,

    pan_dragging: bool,
    last_mouse_screen: Vec2,

    // boom placement drag
    boom_dragging: bool,
    boom_start: Vec2,
}

impl App {
    pub fn new() -> Self {
        let world = World::new(180, 110);

        let mut camera = Camera2D::from_display_rect(Rect::new(0.0, 0.0, 1280.0, 720.0));
        camera.target = vec2(world.w as f32 * world.cell * 0.5, world.h as f32 * world.cell * 0.5);

        Self {
            world,
            controls: Controls {
                paused: false,
                sim_speed: 1.0,

                show_vectors: false,
                show_particles: true,
                show_temp: true,
                show_oil: true,
                show_streamlines: true,
                show_response: true,

                wind_dir_deg: 25.0,
                wind_strength: 1.3,

                advection: 1.0,
                diffusion: 0.04,

                oil_diffusion: 0.08,
                oil_decay: 0.035,
                oil_source_rate: 0.55,
                oil_max: 1.0,

                tool: ToolMode::SpillOil,
                brush_radius: 46.0,

                sl_spacing: 56.0,
                sl_steps: 26,
                sl_dt: 0.016,

                boom_blocks: 0.92,
                skimmer_rate: 1.25,
                dispersant_strength: 1.2,
            },
            camera,
            pan_dragging: false,
            last_mouse_screen: vec2(0.0, 0.0),

            boom_dragging: false,
            boom_start: vec2(0.0, 0.0),
        }
    }

    pub fn handle_camera_input(&mut self) {
        let (_wx, wy) = mouse_wheel();
        if wy.abs() > 0.0 {
            let zoom_mul = (1.0 + (-wy * 0.08)).clamp(0.75, 1.35);
            self.camera.zoom *= zoom_mul;

            self.camera.zoom.x = self.camera.zoom.x.clamp(0.0005, 0.02);
            self.camera.zoom.y = self.camera.zoom.y.clamp(-0.02, -0.0005);
        }

        let m = vec2(mouse_position().0, mouse_position().1);

        if is_mouse_button_pressed(MouseButton::Right) {
            self.pan_dragging = true;
            self.last_mouse_screen = m;
        }
        if !is_mouse_button_down(MouseButton::Right) {
            self.pan_dragging = false;
        }

        if self.pan_dragging {
            let a = self.camera.screen_to_world(self.last_mouse_screen);
            let b = self.camera.screen_to_world(m);
            self.camera.target += a - b;
            self.last_mouse_screen = m;
        }

        let maxx = self.world.w as f32 * self.world.cell;
        let maxy = self.world.h as f32 * self.world.cell;
        self.camera.target.x = self.camera.target.x.clamp(-200.0, maxx + 200.0);
        self.camera.target.y = self.camera.target.y.clamp(-200.0, maxy + 200.0);
    }

    pub fn handle_tool_input(&mut self, ui_wants_pointer: bool) {
        if ui_wants_pointer {
            self.boom_dragging = false;
            return;
        }

        let world_pos = self.camera.screen_to_world(vec2(mouse_position().0, mouse_position().1));

        match self.controls.tool {
            ToolMode::PaintHeat => {
                if is_mouse_button_down(MouseButton::Left) {
                    self.world.splat_temp(world_pos, self.controls.brush_radius, 0.9);
                }
            }
            ToolMode::SpillOil => {
                if is_mouse_button_down(MouseButton::Left) {
                    self.world.splat_oil(world_pos, self.controls.brush_radius, 0.9);
                }
            }
            ToolMode::PlaceSkimmer => {
                if is_mouse_button_pressed(MouseButton::Left) {
                    self.world.skimmers.push(Skimmer { p: world_pos, r: self.controls.brush_radius, rate: self.controls.skimmer_rate });
                }
            }
            ToolMode::PlaceDispersant => {
                if is_mouse_button_pressed(MouseButton::Left) {
                    self.world.dispersants.push(Dispersant { p: world_pos, r: self.controls.brush_radius, strength: self.controls.dispersant_strength });
                }
            }
            ToolMode::EraseResponse => {
                if is_mouse_button_pressed(MouseButton::Left) {
                    self.world.erase_response_near(world_pos, self.controls.brush_radius);
                }
            }
            ToolMode::PlaceBoom => {
                if is_mouse_button_pressed(MouseButton::Left) {
                    self.boom_dragging = true;
                    self.boom_start = world_pos;
                }
                if is_mouse_button_released(MouseButton::Left) && self.boom_dragging {
                    self.boom_dragging = false;
                    let a = self.boom_start;
                    let b = world_pos;
                    if a.distance(b) > 10.0 {
                        self.world.booms.push(BoomSeg { a, b, block: self.controls.boom_blocks });
                    }
                }
            }
        }

        // nice preview while dragging boom
        if self.controls.tool == ToolMode::PlaceBoom && self.boom_dragging {
            set_camera(&self.camera);
            draw_line(self.boom_start.x, self.boom_start.y, world_pos.x, world_pos.y, 3.0, Color::new(1.0, 0.85, 0.2, 0.9));
            set_default_camera();
        }
    }
}

pub fn draw_ui(app: &mut App) -> bool {
    let mut wants_pointer = false;

    egui_macroquad::ui(|ctx| {
        wants_pointer = ctx.wants_pointer_input();

        egui::Window::new("oil spill response")
            .default_pos([12.0, 12.0])
            .resizable(true)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    if ui.button(if app.controls.paused { "resume" } else { "pause" }).clicked() {
                        app.controls.paused = !app.controls.paused;
                    }
                    if ui.button("reset").clicked() {
                        let (w, h) = (app.world.w, app.world.h);
                        app.world = World::new(w, h);
                    }
                });

                ui.add(egui::Slider::new(&mut app.controls.sim_speed, 0.05..=3.0).text("sim speed"));

                ui.separator();

                ui.label("layers");
                ui.checkbox(&mut app.controls.show_temp, "temperature");
                ui.checkbox(&mut app.controls.show_oil, "oil slick");
                ui.checkbox(&mut app.controls.show_particles, "particles");
                ui.checkbox(&mut app.controls.show_streamlines, "streamlines");
                ui.checkbox(&mut app.controls.show_vectors, "vector arrows");
                ui.checkbox(&mut app.controls.show_response, "response assets");

                if app.controls.show_streamlines {
                    ui.add(egui::Slider::new(&mut app.controls.sl_spacing, 24.0..=110.0).text("sl spacing"));
                    ui.add(egui::Slider::new(&mut app.controls.sl_steps, 8..=70).text("sl steps"));
                    ui.add(egui::Slider::new(&mut app.controls.sl_dt, 0.004..=0.04).text("sl dt"));
                }

                ui.separator();

                ui.label("currents");
                ui.add(egui::Slider::new(&mut app.controls.wind_dir_deg, 0.0..=360.0).text("wind dir (deg)"));
                ui.add(egui::Slider::new(&mut app.controls.wind_strength, 0.0..=4.0).text("wind strength"));
                ui.add(egui::Slider::new(&mut app.controls.advection, 0.0..=2.0).text("advection (oil+temp)"));
                ui.add(egui::Slider::new(&mut app.controls.diffusion, 0.0..=0.20).text("temp diffusion"));

                ui.separator();

                ui.label("oil behavior");
                ui.add(egui::Slider::new(&mut app.controls.oil_diffusion, 0.0..=0.35).text("oil diffusion"));
                ui.add(egui::Slider::new(&mut app.controls.oil_decay, 0.0..=0.20).text("weathering/decay"));
                ui.add(egui::Slider::new(&mut app.controls.oil_source_rate, 0.0..=2.0).text("source rate"));
                ui.add(egui::Slider::new(&mut app.controls.oil_max, 0.5..=3.0).text("oil max (visual)"));

                ui.separator();

                ui.label("tools");
                ui.horizontal_wrapped(|ui| {
                    ui.selectable_value(&mut app.controls.tool, ToolMode::SpillOil, "spill oil");
                    ui.selectable_value(&mut app.controls.tool, ToolMode::PlaceBoom, "boom");
                    ui.selectable_value(&mut app.controls.tool, ToolMode::PlaceSkimmer, "skimmer");
                    ui.selectable_value(&mut app.controls.tool, ToolMode::PlaceDispersant, "dispersant");
                    ui.selectable_value(&mut app.controls.tool, ToolMode::EraseResponse, "erase");
                    ui.selectable_value(&mut app.controls.tool, ToolMode::PaintHeat, "paint heat");
                });

                ui.add(egui::Slider::new(&mut app.controls.brush_radius, 8.0..=140.0).text("radius"));

                ui.add(egui::Slider::new(&mut app.controls.boom_blocks, 0.0..=0.98).text("boom block %"));
                ui.add(egui::Slider::new(&mut app.controls.skimmer_rate, 0.2..=3.0).text("skimmer rate"));
                ui.add(egui::Slider::new(&mut app.controls.dispersant_strength, 0.2..=3.0).text("dispersant strength"));

                ui.separator();

                let m = app.world.metrics();
                ui.label(format!("oil mass: {:.2} | max: {:.2} | assets: booms {} skimmers {} dispersants {}",
                                 m.oil_mass, m.oil_max, app.world.booms.len(), app.world.skimmers.len(), app.world.dispersants.len()
                ));

                ui.label("controls: left = use tool | right-drag pan | wheel zoom");
            });
    });

    egui_macroquad::draw();
    wants_pointer
}

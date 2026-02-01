mod render;
mod sim;
mod ui;

use macroquad::prelude::*;

fn conf() -> Conf {
    Conf {
        window_title: "ocean sandbox".to_string(),
        window_width: 1280,
        window_height: 720,
        high_dpi: true,
        ..Default::default()
    }
}

#[macroquad::main(conf)]
async fn main() {
    let mut app = ui::App::new();

    loop {
        clear_background(BLACK);

        app.handle_camera_input();

        if !app.controls.paused {
            let dt = get_frame_time() * app.controls.sim_speed;
            app.world.step(dt, &app.controls);
        }

        render::draw_world(&app.world, &app.controls, &app.camera);

        let ui_wants_pointer = ui::draw_ui(&mut app);
        app.handle_tool_input(ui_wants_pointer);

        next_frame().await;
    }
}

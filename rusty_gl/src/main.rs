use rusty_gl::rusty_gl::{*};
use rusty_gl::sdl2_wrapper::{init_sdl_context, init_video_subsystem, init_window, init_canvas, init_event_pump, Lines};
use sdl2;
use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::rect::Point;
use std::time::Duration;
use rand::prelude::*;

fn main() {
    let mut points: Buffer = Buffer::init();
    let camera = Camera::init(Point3d { x: (0.0), y: (0.0), z: (0.0) }, 1000, 1000, 120.0, 0.1, 5000.0);
    let mut obj = Obj3d::init(
        &[
            Point3d::init(-0.5, -0.5, 0.11),
            Point3d::init(-0.5, 0.5, 0.11),
            Point3d::init(0.5, -0.5, 0.11),
            Point3d::init(0.5, 0.5, 0.11),

            Point3d::init(-0.5, -0.5, 0.12),
            Point3d::init(-0.5, 0.5, 0.12),
            Point3d::init(0.5, -0.5, 0.12),
            Point3d::init(0.5, 0.5, 0.12),
            ],
        &[
                0, 1, 3, 2, 0,
                4, 5, 7, 6, 4,
                6, 2, 3, 7, 5, 1
            ],
            Color::RGB(255, 1, 1),
            Color::RGB(1, 255, 1),
            Color::RGB(1, 1, 255)
    );
    points.load_obj_to_buffer(&obj, &camera);
    println!("{:#?}", points);
    
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem
        .window("rust-sdl2 demo", 1000, 1000)
        .position_centered()
        .build()
        .unwrap();

    let mut canvas = window.into_canvas().build().unwrap();

    let mut event_pump = sdl_context.event_pump().unwrap();
    let mut i = 0;
    let mut rng = rand::thread_rng();
    'running: loop {
        obj.rotate_x_axis(0.1);
        points = Buffer::init();
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit { .. }
                | Event::KeyDown {
                    ..
                } => break 'running,
                _ => {}
            }
        }
        canvas.set_draw_color(Color::RGB(0, 0, 0));
        canvas.clear();
        i = (i + 1) % 255;
        canvas.set_draw_color(Color::RGB(i, 64, 255 - i));
        let (w, h) = canvas.output_size().unwrap();
        
        points.load_obj_to_buffer(&obj, &camera);
        // For performance, it's probably better to draw a whole bunch of points at once
        canvas.draw_lines_w(&points);
        canvas.draw_points_w(&points);
        //canvas.draw_points(points.points.as_slice()).unwrap();
        canvas.present();
        ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60)); // sloppy FPS limit
    }
}

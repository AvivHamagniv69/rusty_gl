use sdl2;
use sdl2::event::Event;
use sdl2::pixels::Color;
use std::time::Duration;

pub mod rusty_gl;
use crate::rusty_gl::Camera;
use crate::rusty_gl::ObjLoader;
use crate::rusty_gl::Point3;
use crate::rusty_gl::Buffer;
use crate::rusty_gl::SdlWrapper;

/*
    TODO:
    fix clipping
    add gui
    add rotation
    add movement
    add camera movement
    add camera rotation
    add triangulation algorithm
 */

fn main() {
    let mut objs: ObjLoader = ObjLoader::init();
    let a = objs.load_obj_file("o.obj");
    println!("{:?}", a);

    let camera = Camera::init(Point3::init(0.0, 0.0, 0.0, None ), 1000, 1000, 70.0, 0.1, 1000.0);
    let mut points: Buffer;

    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem
        .window("rust-sdl2 demo", 1000, 1000)
        .position_centered()
        .build()
        .unwrap();

    let mut canvas: sdl2::render::Canvas<sdl2::video::Window> = window.into_canvas().build().unwrap();

    let mut event_pump = sdl_context.event_pump().unwrap();
    'running: loop {
        //objs.0[0].rotate_y(1.0);

        points = Buffer::init_buffer(&camera);
        points.load_mesh(&objs.0[0]);
        //points.load_mesh(&objs.0[1]);

        for event in event_pump.poll_iter() {
            match event {
                Event::Quit { .. }
                | Event::KeyDown {
                    ..
                } => break 'running,
                _ => {}
            }
        }
        // clear canvas
        canvas.set_draw_color(Color::RGB(0, 0, 0));
        canvas.clear();
                        
        // For performance, it's probably better to draw a whole bunch of points at once
        canvas.draw_all(&points, &camera);
        canvas.draw_triangles(&points, &camera);
        canvas.draw_lines_w(&points, &camera);
        canvas.draw_points_w(&points, &camera);
        
        canvas.present();
        ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60)); // sloppy FPS limit
    }
}

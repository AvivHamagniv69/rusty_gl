use sdl2;
use sdl2::event::Event;
use sdl2::pixels::Color;
use std::time::Duration;

pub mod rusty_gl;
use crate::rusty_gl::*;

/*
    TODO:
    refactor (a lot of things in one object depend on different objects which can create problems in the future when trying to add new features)
    add depth buffer
    add console interface
    add gui
    add camera rotation
    add lighting system
    fix clippings
 */
/*
    axis relative to blender:
    -y is z,

 */
/*
    what cant you do?
    create animations
    create events
    create movement
    create games
 */

fn main() {
    /*
    let mut objs: ObjLoader = ObjLoader {0: Vec::new()};
    objs.load_obj_file("amongus.obj").unwrap();

    let mut camera = Camera::init(Point3::init(0.0, 0.0, 0.0 ), 1000, 1000, 120.0, 1.0, 100.0);
    camera.move_camera(0.0, 0.0, 0.0);
    let mut points: Buffer;

    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem
        .window("rust-sdl2 demo", 1000, 1000).resizable()
        .position_centered()
        .build()
        .unwrap();

    let mut canvas: sdl2::render::Canvas<sdl2::video::Window> = window.into_canvas().build().unwrap();

    let mut event_pump = sdl_context.event_pump().unwrap();
    'running: loop {
        //camera.move_camera(0.0, 0.0, 0.1);
        objs.0[0].rotate_y(0.5);

        points = Buffer {amt_of_points: 0, points: Vec::new(), triangles: Vec::new()};
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
        //canvas.draw_all(&points, &camera);
        canvas.draw_triangles(&mut points, &mut camera).unwrap();
        //canvas.draw_lines_w(&mut points, &camera);
        
        canvas.present();
        ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60)); // sloppy FPS limit
    }
    */
}

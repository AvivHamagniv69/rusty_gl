use sdl2::sys::KeyCode;
use sdl2::{self, keyboard};
use sdl2::event::Event;
use sdl2::pixels::Color;
use std::time::Duration;
use sdl2::keyboard::Scancode;

pub mod rusty_gl;
use crate::rusty_gl::*;

/*
    TODO:
    refactor (a lot of things in one object depend on different objects which can create problems in the future when trying to add new features)
    add depth buffer
    add camera rotation
    add lighting system
    fix clippings
    implement shaders
    add per pixel lighting
    add normal support in obj parser
    add directional lighting
    add multitreading
    fix backface culling errors
 */
/*
    axis relative to blender:
    -y is z,

 */

fn main() {
    let lights = vec![
        Light::Point(PointLight { origin_point: Point3 { x: 1.0, y: 0.0, z: 1.0 }, brightness: 1.0, color: Color { r: 255, g: 0, b: 0, a: 255 }, coefficients: (1.0, 1.0, 1.0) }),
        Light::Point(PointLight { origin_point: Point3 { x: -1.0, y: 0.0, z: 1.0 }, brightness: 1.0, color: Color { r: 255, g: 255, b: 255, a: 255 }, coefficients: (1.0, 1.0, 1.0) }),
        Light::Point(PointLight { origin_point: Point3 { x: 0.0, y: 0.5, z: 1.3 }, brightness: 1.0, color: Color { r: 0, g: 0, b: 255, a: 255 }, coefficients: (1.0, 1.0, 1.0) }),
    ];
    let mut objs: MeshLoader = MeshLoader {0: Vec::new()};
    objs.load_obj_file("sphere.obj").unwrap();

    let mut camera = Camera::init(Point3 {x: 0.0, y: 0.0, z: 0.0 }, 1000, 1000, 80.0, 1.0, 100.0);
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
        objs.0[0].set_origin_point();
        objs.0[0].rotate_y(0.5);

        points = Buffer::init();
        points.load_meshes(&objs.0);

        for event in event_pump.poll_iter() {
            match event {
                Event::Quit { .. } => break 'running,
                
                Event::KeyDown { timestamp, window_id, keycode, scancode, keymod, repeat } => {
                    if keycode == Some(keyboard::Keycode::A) {
                        camera.origin_point.x += 0.08;
                    }
                    else if keycode == Some(keyboard::Keycode::D) {
                        camera.origin_point.x -= 0.08;
                    }
                    else if keycode == Some(keyboard::Keycode::W) {
                        camera.origin_point.z += 0.08;
                    }
                    else if keycode == Some(keyboard::Keycode::S) {
                        camera.origin_point.z -= 0.08;
                    }
                    else if keycode == Some(keyboard::Keycode::Space) {
                        camera.origin_point.y += 0.08;
                    }
                    else if keycode == Some(keyboard::Keycode::CapsLock) {
                        camera.origin_point.y -= 0.08;
                    }
                    else if keycode == Some(keyboard::Keycode::E) {
                        camera.rotate_y(1.5);
                    }
                    else if keycode == Some(keyboard::Keycode::Q) {
                        camera.rotate_y(-1.5);
                    }
                    else if keycode == Some(keyboard::Keycode::F) {
                        camera.rotate_x(1.5);
                    }
                    else if keycode == Some(keyboard::Keycode::C) {
                        camera.rotate_x(-1.5);
                    }
                }
                _ => {}
            }
        }

        // clear canvas
        canvas.set_draw_color(Color::RGB(0, 0, 0));
        canvas.clear();
                        
        // For performance, it's probably better to draw a whole bunch of points at once
        //canvas.draw_all(&points, &camera);
        canvas.draw_triangles(&mut points, &mut camera, &lights).unwrap();
        //canvas.draw_lines_w(&mut points, &camera);
        
        canvas.present();
        ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60)); // sloppy FPS limit
    }
}

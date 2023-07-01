use sdl2::event::Event;
use sdl2::pixels::Color;
use sdl2::render::Canvas;
use sdl2::video::Window;
use sdl2::{self, keyboard};
pub mod rusty_gl;
use crate::rusty_gl::*;

/*
   TODO:
   add depth buffer
   add camera rotation
   add lighting system
   fix clippings
   implement shaders
   add per pixel lighting
   add normal support in obj parser
   add directional lighting
   add multithreading
   add async
   fix backface culling errors
*/
/*
   axis relative to blender:
   -y is z,

*/

fn main() {
    RUSTY_GL.set_amt_of_threads(15);
    let lights = vec![Light::Point(PointLight {
        origin_point: Point3 {
            x: 1.0,
            y: 0.0,
            z: 1.3,
        },
        brightness: 10.0,
        color: Color {
            r: 255,
            g: 0,
            b: 0,
            a: 255,
        },
        coefficients: (1.0, 1.0, 1.0),
    })];
    let mut objs: MeshLoader = MeshLoader { 0: Vec::new() };
    objs.load_obj_file("cube.obj").unwrap();

    let mut camera = Camera::init(
        Point3 {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        },
        1000,
        1000,
        80.0,
        1.0,
        100.0,
    );
    let mut points: Buffer;

    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window: Window = video_subsystem
        .window("rusty_gl", 1000, 1000)
        .resizable()
        .position_centered()
        .build()
        .unwrap();

    let mut canvas = CanvasWrapper::init(window);

    let mut event_pump = sdl_context.event_pump().unwrap();
    'running: loop {
        objs.0[0].set_origin_point();
        objs.0[0].rotate_y(0.05);

        points = Buffer::init();
        points.load_meshes(&objs.0);

        for event in event_pump.poll_iter() {
            match event {
                Event::Quit { .. } => break 'running,

                Event::KeyDown {
                    timestamp,
                    window_id,
                    keycode,
                    scancode,
                    keymod,
                    repeat,
                } => {
                    if keycode == Some(keyboard::Keycode::A) {
                        camera.origin_point.x += 0.08;
                    } else if keycode == Some(keyboard::Keycode::D) {
                        camera.origin_point.x -= 0.08;
                    } else if keycode == Some(keyboard::Keycode::W) {
                        camera.origin_point.z += 0.08;
                    } else if keycode == Some(keyboard::Keycode::S) {
                        camera.origin_point.z -= 0.08;
                    } else if keycode == Some(keyboard::Keycode::Space) {
                        camera.origin_point.y += 0.08;
                    } else if keycode == Some(keyboard::Keycode::CapsLock) {
                        camera.origin_point.y -= 0.08;
                    } else if keycode == Some(keyboard::Keycode::E) {
                        camera.rotate_y(1.5);
                    } else if keycode == Some(keyboard::Keycode::Q) {
                        camera.rotate_y(-1.5);
                    } else if keycode == Some(keyboard::Keycode::F) {
                        camera.rotate_x(1.5);
                    } else if keycode == Some(keyboard::Keycode::C) {
                        camera.rotate_x(-1.5);
                    }
                }
                _ => {}
            }
        }

        canvas.0.set_draw_color(Color::RGB(0, 0, 0));
        canvas.0.clear();

        canvas.draw_triangles(points, &mut camera, &lights).unwrap();
        canvas.0.present();
    }
}

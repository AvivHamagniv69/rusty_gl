use sdl2;
use sdl2::render::Canvas;
use sdl2::video::Window;

use crate::rusty_gl::{Buffer};

pub fn init_sdl_context() -> sdl2::Sdl {
    return sdl2::init().unwrap();
}
pub fn init_video_subsystem(sdl_context: sdl2::Sdl) -> sdl2::VideoSubsystem {
    return sdl_context.video().unwrap();
}
pub fn init_window(video_subsystem: sdl2::VideoSubsystem, title: &str, width: u32, height: u32) -> sdl2::video::Window {
    return video_subsystem.window(title, width, height).build().unwrap()
}
pub fn init_canvas(window: sdl2::video::Window) -> sdl2::render::Canvas<sdl2::video::Window> {
    return window.into_canvas().build().unwrap();
}
pub fn init_event_pump(sdl_context: sdl2::Sdl) -> sdl2::EventPump {
    return sdl_context.event_pump().unwrap();
}

pub trait Lines {
    fn draw_lines_w(&mut self, buffer: &Buffer);
    fn draw_points_w(&mut self, buffer: &Buffer);
}

impl Lines for Canvas<Window> {
    fn draw_lines_w(&mut self, buffer: &Buffer) {
        for i in 0..(*buffer).indices.len()-1 {
            self.set_draw_color(buffer.lines_colors[i]);
            self.draw_line(buffer.points[buffer.indices[i]], buffer.points[buffer.indices[i+1]]).unwrap();
        }
    }
    fn draw_points_w(&mut self, buffer: &Buffer) {
        for i in 0..buffer.points.len() {
            self.set_draw_color(buffer.points_colors[i]);
            self.draw_point(buffer.points[i]).unwrap();
        }
    }
}
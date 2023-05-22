use sdl2::{self, rect::Point};
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

pub trait Triangles {
    fn draw_triangle(&mut self, points: [Point; 3]);
    fn draw_triangles(&mut self, buffer: &Buffer);
}

impl Triangles for Canvas<Window> {
    fn draw_triangle(&mut self, points: [Point; 3]) {
        let p_min = 
        if points[0].y <= points[1].y && points[0].y <= points[2].y {
            points[0]
        }
        else if points[1].y <= points[0].y && points[1].y <= points[2].y {
            points[1]
        }
        else {
            points[2]
        };
        
        let p_max = 
        if points[0].y >= points[1].y && points[0].y >= points[2].y {
            points[0]
        }
        else if points[1].y >= points[0].y && points[1].y >= points[2].y {
            points[1]
        } else {
            points[2]
        };

        let p_mid = 
        if points[0] != p_max && points[0] != p_min {
            points[0]
        }
        else if points[1] != p_max && points[1] != p_min {
            points[1]
        }
        else {
            points[2]
        };

        let mut m_min_to_max: f64;
        if p_min.x as f64 - p_max.y as f64 != 0.0 && p_min.x as f64 - p_max.y as f64 != -0.0 {
            m_min_to_max = (p_min.y as f64 - p_max.y as f64) / (p_min.x as f64 - p_max.y as f64);
        }  
        else if p_min.x as f64 - p_max.y as f64 == -0.0 {
            m_min_to_max = (p_min.y as f64 - p_max.y as f64) / (-0.00000000001 as f64);
        }
        else {
            m_min_to_max = (p_min.y as f64 - p_max.y as f64) / (0.00000000001 as f64);
        }
        let b_min_to_max: f64 = -m_min_to_max * p_min.x as f64 + p_min.y as f64;

        let mut m_min_to_mid: f64;
        if p_min.x as f64 - p_mid.y as f64 != 0.0 && p_min.x as f64 - p_mid.y as f64 != -0.0 {
            m_min_to_mid = (p_min.y as f64 - p_mid.y as f64) / (p_min.x as f64 - p_mid.y as f64);
        }
        else if p_min.x as f64 - p_mid.y as f64 == -0.0 {
            m_min_to_mid = (p_min.y as f64 - p_mid.y as f64) / (-0.00000000001 as f64);
        }
        else {
            m_min_to_mid = (p_min.y as f64 - p_mid.y as f64) / (0.00000000001 as f64);
        }
        let b_min_to_mid: f64 = -m_min_to_mid * p_min.x as f64 + p_min.y as f64;

        let mut m_mid_to_max: f64;
        if p_max.x as f64 - p_mid.y as f64 != 0.0 && p_max.x as f64 - p_mid.y as f64 != -0.0 {
            m_mid_to_max = (p_max.y as f64 - p_mid.y as f64) / (p_max.x as f64 - p_mid.y as f64);
        }  
        else if p_max.x as f64 - p_mid.y as f64 == -0.0 {
            m_mid_to_max = (p_max.y as f64 - p_mid.y as f64) / (-0.00000000001 as f64)
        }
        else {
            m_mid_to_max = (p_max.y as f64 - p_mid.y as f64) / (0.00000000001 as f64)
        }
        
        let b_mid_to_max: f64 = -m_mid_to_max * p_mid.x as f64 + p_mid.y as f64;

        println!("{:#?}", p_min);
        println!("{:#?}", p_mid);

        for y in p_min.y..p_mid.y {
            let x_start = (y as f64/m_min_to_max) - (b_min_to_max/m_min_to_max);
            let x_end = (y as f64/m_min_to_mid) - (b_min_to_mid/m_min_to_mid);
            self.draw_line(Point::new(x_start as i32, y), Point::new(x_end as i32, y)).unwrap();
        }
        for y in p_mid.y..p_max.y {
            let x_start = (y as f64/m_min_to_max) - (b_min_to_max/m_min_to_max);
            let x_end = (y as f64/m_mid_to_max) - (b_mid_to_max/m_mid_to_max);
            self.draw_line(Point::new(x_start as i32, y), Point::new(x_end as i32, y)).unwrap();
        }
    }

    fn draw_triangles(&mut self, buffer: &Buffer) {
        
    }
}

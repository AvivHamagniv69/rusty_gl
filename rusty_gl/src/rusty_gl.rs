use sdl2::render::Canvas;
use std::f32::consts::PI;
use std::str;
use std::fs::File;
use std::error::Error;
use std::io::BufReader;
use std::io::BufRead;
use std::sync::Arc;
use std::time::Duration;

use sdl2;
use sdl2::pixels::Color;
use sdl2::rect::Point;
use sdl2::video::Window;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Point3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}
impl Point3 {
    pub fn init(x: f32, y: f32, z: f32) -> Point3 {
        Point3 { x: (x), y: (y), z: (z)}
    }    

    fn project_point(point: Self, camera: &Camera) -> Option<Point3> {
        // a = forward
        // b = right
        // c = up
        let perspective = [
            1.0/((camera.width as f32/camera.height as f32)*f32::tan(camera.fov_radians/2.0)),
            1.0/f32::tan(camera.fov_radians/2.0),
            -((camera.far_clipping_plane+camera.near_clipping_plane)/(camera.far_clipping_plane-camera.near_clipping_plane)),
            -((2.0*camera.far_clipping_plane*camera.near_clipping_plane)/(camera.far_clipping_plane-camera.near_clipping_plane))
        ];

        let x: f32 = point.x-camera.origin_point.x;
        let y: f32 = point.y-camera.origin_point.y;
        let z: f32 = point.z-camera.origin_point.z;
        let w: f32 = -z;
        
        //let temp_x = x;
        //let temp_y = y;
        //let temp_z = z;
        //let temp_w = w;

        // checks if the point is behind the camera
        if z < camera.near_clipping_plane || z > camera.far_clipping_plane {
            return None;
        }
        return Some(Point3::init(
            x * perspective[0] / w * camera.fov_radians * camera.width as f32 + camera.width as f32/2.0,
            y * perspective[1] / w * camera.fov_radians * camera.height as f32 + camera.height as f32/2.0,
            w
        ));
    }

    pub fn calculate_distance(point1: &Self, point2: &Self) -> f32 {
        f32::sqrt((point1.x-point2.x).powi(2) + (point1.y-point2.y).powi(2) + (point1.z-point2.z).powi(2))
    }

    pub fn get_dot_product(point1: &Point3, point2: &Point3) -> f32 {
        point1.x * point2.x + point1.y * point2.y + point1.z + point2.z
    }

    pub fn normalize_vector(point: Point3) -> Point3 {
        let m = (point.x.powi(2) + point.y.powi(2) + point.z.powi(2)).sqrt();

        Point3 { 
            x: (point.x/m),
            y: (point.y/m),
            z: (point.z/m)
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Triangle {
    pub point1: usize,
    pub point2: usize,
    pub point3: usize,

    origin_point: Option<Point3>,
    normal: Option<Point3>,

    pub color: Option<Color>,
}

impl Triangle {
    pub fn calc_origin_point_triangle(points: [Point3; 3]) -> Point3 {
        Point3 { 
            x: (points[0].x+points[1].x+points[2].x)/3.0,
            y: (points[0].y+points[1].y+points[2].y)/3.0,
            z: (points[0].z+points[1].z+points[2].z)/3.0
        }
    }
    
    pub fn calc_normal(points: [Point3; 3]) -> Point3 {
        let vu = Point3::init(points[1].x-points[0].x, points[1].y-points[0].y, points[1].z-points[0].z);
        let vv = Point3::init(points[2].x-points[0].x, points[2].y-points[0].y, points[2].z-points[0].z);
    
        let a = Point3 {
            x: (vu.y*vv.z) - (vu.z*vv.y),
            y: (vu.z*vv.x) - (vu.x*vv.z),
            z: (vu.x*vv.y) - (vu.y*vv.x)
        };
        Point3::normalize_vector(a)
    }
}

pub trait MeshAndBuffer {
    fn sort_triangles_in_mesh_origin_point_z_distance_from_camera(&mut self, camera: &Camera) -> Result<(), ()>;
}

#[derive(Debug, Clone)]
pub struct Mesh {
    pub name: String,

    pub points: Vec<Point3>,
    pub triangles: Vec<Triangle>,

    pub origin_point: Point3,
    pub color: Option<Color>
}

impl MeshAndBuffer for Mesh {
    fn sort_triangles_in_mesh_origin_point_z_distance_from_camera(&mut self, camera: &Camera) -> Result<(), ()> {
        if self.triangles.len() == 0 || self.triangles.len() == 1 {return Err(());}
        
        for i in 0..self.triangles.len() {
            self.triangles[i].origin_point.unwrap().x -= camera.origin_point.x;
            self.triangles[i].origin_point.unwrap().y -= camera.origin_point.y;
            self.triangles[i].origin_point.unwrap().z -= camera.origin_point.z;
        }

        for i in 1..self.triangles.len() {
            if self.triangles[i-1].origin_point.unwrap().z > self.triangles[i].origin_point.unwrap().z {
                self.triangles.sort_by(|a, b| a.origin_point.unwrap().z.partial_cmp(&b.origin_point.unwrap().z).unwrap());
                break;
            }
        }
        
        for i in 0..self.triangles.len() {
            self.triangles[i].origin_point.unwrap().x += camera.origin_point.x;
            self.triangles[i].origin_point.unwrap().x += camera.origin_point.y;
            self.triangles[i].origin_point.unwrap().x += camera.origin_point.z;
        }
        Ok(())
    }
}

impl Mesh {
    pub fn init(name: String, points: &[Point3], triangles: &[Triangle], color: Color) -> Mesh {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        for i in points {
            x += i.x;
            y += i.y;
            z += i.z;
        }
        let p = Point3 {x: (x/points.len() as f32), y: (y/points.len() as f32), z: (z/points.len() as f32)};
        Mesh {name: (name), points: (points.to_vec()), triangles: (triangles.to_vec()), origin_point: (p), color: (Some(color))}
    }

    pub fn set_origin_point(&mut self) {
        let mut p = Point3::init(0.0, 0.0, 0.0);
        for i in 0..self.points.len() {
            p.x += self.points[i].x;
            p.y += self.points[i].y;
            p.z += self.points[i].z;
        }
        p.x /= self.points.len() as f32;
        p.y /= self.points.len() as f32;
        p.z /= self.points.len() as f32;
        self.origin_point = p;
    }

    pub fn rotate_x(&mut self, mut angle: f32) {
        angle *= PI/180.0;
        self.set_origin_point();
        for p in 0..self.points.len() {
            self.points[p].y -= self.origin_point.y;
            self.points[p].z -= self.origin_point.z;

            let y = self.points[p].y;
            let z = self.points[p].z;

            self.points[p].y = y * f32::cos(angle) - z * f32::sin(angle);
            self.points[p].z = y * f32:: sin(angle) + z * f32::cos(angle);

            self.points[p].y += self.origin_point.y;
            self.points[p].z += self.origin_point.z;
        }
        self.set_origin_point();
        self.set_normals();
    }

    pub fn rotate_y(&mut self, mut angle: f32) {
        angle *= PI/180.0;
        self.set_origin_point();
        for p in 0..self.points.len() {
            self.points[p].x -= self.origin_point.x;
            self.points[p].z -= self.origin_point.z;

            let x = self.points[p].x;
            let z = self.points[p].z;

            self.points[p].x = x * f32::cos(angle) + z * f32::sin(angle);
            self.points[p].z = -(x * f32:: sin(angle)) + z * f32::cos(angle);

            self.points[p].x += self.origin_point.x;
            self.points[p].z += self.origin_point.z;
        }
        self.set_origin_point();
        self.set_normals();
    }

    pub fn rotate_z(&mut self, mut angle: f32) {
        angle *= PI/180.0;

        self.set_origin_point();
        for p in 0..self.points.len() {
            self.points[p].x -= self.origin_point.x;
            self.points[p].y -= self.origin_point.y;

            let x = self.points[p].x;
            let y = self.points[p].y;

            self.points[p].x = x * f32::cos(angle) - y * f32::sin(angle);
            self.points[p].y = x * f32:: sin(angle) + y * f32::cos(angle);

            self.points[p].x += self.origin_point.x;
            self.points[p].y += self.origin_point.y;
        }
        self.set_origin_point();
        self.set_normals();
    }

    pub fn set_normals(&mut self) {
        for i in 0..self.triangles.len() {
            self.triangles[i].normal = Some(Triangle::calc_normal([self.points[self.triangles[i].point1], self.points[self.triangles[i].point2], self.points[self.triangles[i].point3]]));
        }
        /*
        for i in 0..self.polygons.len() {
            let mut n = Point3 {x: 0.0, y: 0.0, z: 0.0};
            for j in 0..self.polygons[i].points.len()-1 {
                n.x += self.points[self.polygons[i].points[j]].x * self.points[self.polygons[i].points[j+1]].x + self.points[self.polygons[i].points[j]].x * self.points[self.polygons[i].points[j+1]].y + self.points[self.polygons[i].points[j]].x * self.points[self.polygons[i].points[j+1]].z;
                n.y += self.points[self.polygons[i].points[j]].y * self.points[self.polygons[i].points[j+1]].x + self.points[self.polygons[i].points[j]].y * self.points[self.polygons[i].points[j+1]].y + self.points[self.polygons[i].points[j]].y * self.points[self.polygons[i].points[j+1]].z;
                n.z += self.points[self.polygons[i].points[j]].z * self.points[self.polygons[i].points[j+1]].x + self.points[self.polygons[i].points[j]].z * self.points[self.polygons[i].points[j+1]].y + self.points[self.polygons[i].points[j]].z * self.points[self.polygons[i].points[j+1]].z;
            }
        }
        */
    }
}

#[derive(Debug, Clone)]
pub struct Buffer {
    pub amt_of_points: usize,

    pub points: Vec<Point3>,
    pub triangles: Vec<Triangle>,
}
impl Buffer {
    pub fn load_mesh(&mut self, mesh: &Mesh) {
        for i in 0..mesh.points.len() {
            self.points.push(mesh.points[i]);
        }
        for i in 0..mesh.triangles.len() {
            let l = Triangle::calc_normal([mesh.points[mesh.triangles[i].point1], mesh.points[mesh.triangles[i].point2], mesh.points[mesh.triangles[i].point3]]);
            self.triangles.push(Triangle { point1: mesh.triangles[i].point1+self.amt_of_points, point2: mesh.triangles[i].point2+self.amt_of_points, point3: mesh.triangles[i].point3+self.amt_of_points, origin_point: mesh.triangles[i].origin_point, color: mesh.triangles[i].color, normal: Some(l) });
        }
        self.amt_of_points += mesh.points.len();
    }
}

impl<'a> MeshAndBuffer for Buffer {
    fn sort_triangles_in_mesh_origin_point_z_distance_from_camera(&mut self, camera: &Camera) -> Result<(), ()> {
        if self.triangles.len() == 0 {return Err(());}
        
        for i in 0..self.triangles.len() {
            self.triangles[i].origin_point.unwrap().x -= camera.origin_point.x;
            self.triangles[i].origin_point.unwrap().y -= camera.origin_point.y;
            self.triangles[i].origin_point.unwrap().z -= camera.origin_point.z;
        }

        self.triangles.sort_by(|a, b| b.origin_point.unwrap().z.partial_cmp(&a.origin_point.unwrap().z).unwrap());
        
        for i in 0..self.triangles.len() {
            self.triangles[i].origin_point.unwrap().x += camera.origin_point.x;
            self.triangles[i].origin_point.unwrap().x += camera.origin_point.y;
            self.triangles[i].origin_point.unwrap().x += camera.origin_point.z;
        }
        Ok(())
    }
}
pub trait Lighting {
    fn apply_basic_shading(color: Color, normal: Point3, light: Point3, coefficient: f32) -> Color {
        let c = Point3::get_dot_product(&light, &normal);
        Color { 
            r: (color.r as f32 * c * coefficient).abs() as u8,
            g: (color.g as f32 * c * coefficient).abs() as u8,
            b: (color.b as f32 * c * coefficient).abs() as u8,
            a: color.a
        }
    }
}

pub struct Light {
    pub origin_point: Point3
}

#[derive(Debug, Clone)]
pub struct Camera {
    pub origin_point: Point3,

    pub width: u32,
    pub height: u32,

    fov: f32,
    fov_radians: f32,

    near_clipping_plane: f32,
    far_clipping_plane: f32,

    z_depth_buffer: Vec<f32>,

    pub forward_vec: Point3,
    pub right_vec: Point3,
    pub up_vec: Point3,
    pub look_vec: Point3,
    pub target_vec: Point3,
}

impl<'a> Lighting for Camera {}

impl<'a> Camera {
    pub fn init(origin_point: Point3, width: u32, height: u32, fov: f32, near_clipping_plane: f32, far_clipping_plane: f32) -> Camera {
        Camera { 
            origin_point: (origin_point),
            width: (width),
            height: (height),
            fov: (fov),
            fov_radians: (fov*PI as f32/180.0),
            near_clipping_plane: (near_clipping_plane),
            far_clipping_plane: (far_clipping_plane),
            z_depth_buffer: vec![far_clipping_plane; width as usize * height as usize],
            forward_vec: Point3::normalize_vector(Point3 { x: 0.0, y: 0.0, z: 1.0 }),
            right_vec: Point3::normalize_vector(Point3 { x: -1.0, y: 0.0, z: 0.0 }),
            up_vec: Point3::normalize_vector(Point3 { x: 0.0, y: 1.0, z: 0.0 }),
            look_vec: Point3::normalize_vector(Point3 { x: 0.0, y: 0.0, z: 1.0 }),
            target_vec: Point3::normalize_vector(Point3 { x: 0.0, y: 0.0, z: 1.0 }),
        }
    }

    pub fn move_camera(&mut self, x: f32, y: f32, z: f32) {
        self.origin_point.x += x;
        self.origin_point.y += y;
        self.origin_point.z += z;
    }
}

pub struct ObjLoader(pub Vec<Mesh>);
impl ObjLoader {
    pub fn load_obj_file(&mut self, path:&str) -> Result<(), Box<dyn Error>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut amt_of_points = 0;
        let mut counter = 0;

        for line in reader.lines() {
            let line_string = line.unwrap();
            let mut line_split = line_string.split(" ");
            let line_split_len = line_split.clone().count();
            let line_0 = line_split.nth(0).unwrap().to_string(); 
            let obj_len = self.0.len();

            if line_0 == "o" {
                if counter > 0 {
                    amt_of_points += self.0[obj_len-1].points.len();
                }
                self.0.push(Mesh::init(line_split.nth(0).unwrap().to_string(), &[],  &[], Color::RGB(180, 180, 180)));
                counter += 1;
            }
            
            if line_0 == "v" {
                let x = line_split.nth(0).unwrap().to_string();
                let y = line_split.nth(0).unwrap().to_string();
                let z = line_split.nth(0).unwrap().to_string();
                self.0[obj_len-1].points.push(Point3::init(x.parse()?, y.parse()?, z.parse()?));
            }

            if line_0 == "f" {
                if line_split_len == 4 {
                    let mut p_1: u32 = line_split.nth(0).unwrap().to_string().parse()?;
                    let mut p_2: u32 = line_split.nth(0).unwrap().to_string().parse()?;
                    let mut p_3: u32 = line_split.nth(0).unwrap().to_string().parse()?;

                    p_1 -= 1 + amt_of_points as u32;
                    p_2 -= 1 + amt_of_points as u32;
                    p_3 -= 1 + amt_of_points as u32;

                    let p = Point3::init(
                        (self.0[obj_len-1].points[p_1 as usize].x + self.0[obj_len-1].points[p_2 as usize].x + self.0[obj_len-1].points[p_2 as usize].x)/3.0,
                        (self.0[obj_len-1].points[p_1 as usize].y + self.0[obj_len-1].points[p_2 as usize].y + self.0[obj_len-1].points[p_2 as usize].y)/3.0,
                        (self.0[obj_len-1].points[p_1 as usize].z + self.0[obj_len-1].points[p_2 as usize].z + self.0[obj_len-1].points[p_2 as usize].z)/3.0
                    );
                    let l = self.0[obj_len-1].color;
                    let n = Triangle::calc_normal([self.0[obj_len-1].points[p_1 as usize], self.0[obj_len-1].points[p_2 as usize], self.0[obj_len-1].points[p_3 as usize]]);
                    self.0[obj_len-1].triangles.push(Triangle { point1: (p_1 as usize), point2: (p_2 as usize), point3: (p_3 as usize), origin_point: Some(p), color: (l), normal: Some(n) });
                }
            }
        }  
        Ok(())
    }
}

pub trait SdlWrapper {
    fn draw_lines_w(&mut self, buffer: &mut Buffer, camera: &Camera) -> Result<(), ()>;
    
    fn draw_triangle(&mut self, points: [Point3; 3], camera: &mut Camera);
    fn draw_triangles(&mut self, buffer: &mut Buffer, camera: &mut Camera) -> Result<(), ()>;

    fn draw_all(&mut self, buffer: &mut Buffer, camera: &Camera);

    fn check_backface_culling(points: [Point3; 3], normal: Point3, camera: &Camera) -> bool {
        normal.x * (points[0].x - camera.origin_point.x) +
        normal.y * (points[0].y - camera.origin_point.y) +
        normal.z * (points[0].z - camera.origin_point.z) > 0.0
    }
}

impl SdlWrapper for Canvas<Window> {
    fn draw_lines_w(&mut self, buffer: &mut Buffer, camera: &Camera) -> Result<(), ()> {
        let check = buffer.sort_triangles_in_mesh_origin_point_z_distance_from_camera(camera);
        if check == Err(()) {return Err(());}
        self.set_draw_color(Color::RGB(255, 255, 255));

        for i in 0..buffer.triangles.len() {
            if buffer.triangles[i].normal == None {return Err(());}

            if Canvas::check_backface_culling([buffer.points[buffer.triangles[i].point1], buffer.points[buffer.triangles[i].point2], buffer.points[buffer.triangles[i].point3]], buffer.triangles[i].normal.unwrap(), camera)
            {
                continue;
            }            
            let p1 = Point3::project_point(buffer.points[buffer.triangles[i].point1], camera);
            let p2 = Point3::project_point(buffer.points[buffer.triangles[i].point2], camera);
            let p3 =  Point3::project_point(buffer.points[buffer.triangles[i].point3], camera);

            if p1 == None || p2 == None || p3 == None {
                continue;
            }

            self.draw_line(Point::new(p1.unwrap().x as i32, p1.unwrap().y as i32), Point::new(p2.unwrap().x as i32, p2.unwrap().y as i32)).unwrap();
            self.draw_line(Point::new(p1.unwrap().x as i32, p1.unwrap().y as i32), Point::new(p3.unwrap().x as i32, p3.unwrap().y as i32)).unwrap();
            self.draw_line(Point::new(p2.unwrap().x as i32, p2.unwrap().y as i32), Point::new(p3.unwrap().x as i32, p3.unwrap().y as i32)).unwrap();
        }
        Ok(())
    }

    fn draw_triangle(&mut self, points_: [Point3; 3], camera: &mut Camera) {
        let mut points = points_;

        if points[1].y < points[0].y {points.swap(0, 1)};
        if points[2].y < points[0].y {points.swap(2, 0)};
        if points[2].y < points[1].y {points.swap(2, 1)};

        let m_0_2;
        let m_0_1;
        let m_1_2;

        if points[2].x-points[0].x != 0.0 {
            m_0_2 = (points[2].y-points[0].y)/(points[2].x-points[0].x);
        }
        else {
            m_0_2 = (points[2].y-points[0].y)/(0.1);
        }

        if points[1].x-points[0].x != 0.0 {
            m_0_1 = (points[1].y-points[0].y)/(points[1].x-points[0].x);
        }
        else {
            m_0_1 = (points[1].y-points[0].y)/(0.1);
        }

        if points[2].x-points[1].x != 0.0 {
            m_1_2 = (points[2].y-points[1].y)/(points[2].x-points[1].x);
        }
        else {
            m_1_2 = (points[2].y-points[1].y)/(0.1);
        }

        let b_0_2 = points[0].y - m_0_2 * points[0].x;
        let b_0_1 = points[0].y - m_0_1 * points[0].x;
        let b_1_2 = points[1].y - m_1_2 * points[1].x;

        for y in points[0].y as i32+1..=points[2].y as i32 {
            if y <= points[1].y as i32 {
                let left_x = -((b_0_1-y as f32)/m_0_1);
                let right_x = -((b_0_2-y as f32)/m_0_2);
                self.draw_line(Point::new(left_x as i32, y), Point::new(right_x as i32, y)).unwrap(); 
            }
            else {
                let left_x = -((b_1_2-y as f32)/m_1_2);
                let right_x = -((b_0_2-y as f32)/m_0_2);
                self.draw_line(Point::new(left_x as i32, y), Point::new(right_x as i32, y)).unwrap(); 
            }
        }
    }

    fn draw_triangles(&mut self, buffer: &mut Buffer, camera: &mut Camera) -> Result<(), ()> {
        //let check = buffer.sort_triangles_in_mesh_origin_point_z_distance_from_camera(camera);
        //if check == Err(()) {return Err(());}

        for t in 0..buffer.triangles.len() {
            if buffer.triangles[t].normal == None {
                return Err(());
            }

            if Canvas::check_backface_culling([buffer.points[buffer.triangles[t].point1], buffer.points[buffer.triangles[t].point2], buffer.points[buffer.triangles[t].point3]], buffer.triangles[t].normal.unwrap(), camera) {
                continue;
            }
            
            if buffer.triangles[t].color == None {
                let color = Camera::apply_basic_shading(Color::RGB(180, 180, 180), buffer.triangles[t].normal.unwrap(), camera.origin_point, 1.0);
                self.set_draw_color(color);
            }
            else {
                let color = Camera::apply_basic_shading(buffer.triangles[t].color.unwrap(), buffer.triangles[t].normal.unwrap(), camera.origin_point, 1.0);
                self.set_draw_color(color);
                
            }
            let p1 = Point3::project_point(buffer.points[buffer.triangles[t].point1], camera);
            let p2 = Point3::project_point(buffer.points[buffer.triangles[t].point2], camera);
            let p3 = Point3::project_point(buffer.points[buffer.triangles[t].point3], camera);

            if p1 == None || p2 == None || p3 == None {
                continue;
            }

            self.draw_triangle([p1.unwrap(), p2.unwrap(), p3.unwrap()], camera);
        }
        camera.z_depth_buffer = vec![camera.far_clipping_plane; camera.width as usize * camera.height as usize];
        Ok(())
    }

    fn draw_all(&mut self, buffer: &mut Buffer, camera: &Camera) {

    }
}

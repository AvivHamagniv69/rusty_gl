use sdl2::render::Canvas;
use std::collections::btree_map::Iter;
use std::f32::consts::PI;
use std::str;
use std::fs::File;
use std::error::Error;
use std::io::BufReader;
use std::io::BufRead;

use sdl2;
use sdl2::pixels::Color;
use sdl2::rect::Point;
use sdl2::video::Window;

pub enum LightingMode {
    PerPixel,
    PerTriangle,
}

pub struct RustyGl {
    lighting_mode: LightingMode,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Point3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}
impl Point3 {
    pub fn rotate_x_around_point_degrees(point: &mut Point3, origin_point: &Point3, angle: f32) {
        let a =angle*PI/180.0;
        
        point.y -= origin_point.y;
        point.z -= origin_point.z;

        let y = point.y;
        let z = point.z;

        point.y = y * a.cos() - z*a.sin();
        point.z = z * a.sin() + z * a.cos();

        point.y += origin_point.y;
        point.z += origin_point.z;
            }

    pub fn rotate_y_around_point_degrees(point: &mut Point3, origin_point: &Point3, angle: f32) {
        let a = angle*PI/180.0;

        point.x -= origin_point.x;
        point.z -= origin_point.z;

        let x = point.x;
        let z = point.z;

        point.x = x * a.cos() + z * a.sin();
        point.z = -(x * a.sin()) + z * a.cos();

        point.x += origin_point.x;
        point.z += origin_point.z;
    }

    pub fn rotate_z_around_point_degrees(point: &mut Point3, origin_point: &Point3, angle: f32) {
        let a = angle*PI/180.0;

        point.x -= origin_point.x;
        point.y -= origin_point.y;

        let x = point.x;
        let y = point.y;

        point.x = x * a.cos() - y * a.sin();
        point.y = x * a.sin() + y * a.cos();

        point.x -= origin_point.x;
        point.y -= origin_point.y;
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

        let mut x: f32 = point.x;
        let mut y: f32 = point.y;
        let mut z: f32 = point.z;
        
        //let temp_x = x;
        //let temp_y = y;
        //let temp_z = z;
        //let temp_w = w;

        let mut temp = Point3 { x: x, y: y, z: z };
        Point3::rotate_x_around_point_degrees(&mut temp, &camera.origin_point, camera.rotations.0);
        Point3::rotate_y_around_point_degrees(&mut temp, &camera.origin_point, camera.rotations.1);
        Point3::rotate_z_around_point_degrees(&mut temp, &camera.origin_point, camera.rotations.2);

        x = temp.x-camera.origin_point.x;
        y = temp.y-camera.origin_point.y;
        z = temp.z-camera.origin_point.z;
        let w = -z;

        if z < camera.near_clipping_plane || z > camera.far_clipping_plane {
            return None;
        }

        let a = Point3 {
            x: x * perspective[0] / w * camera.fov_radians * camera.width as f32 + camera.width as f32/2.0,
            y: y * perspective[1] / w * camera.fov_radians * camera.height as f32 + camera.height as f32/2.0,
            z: w
        };

        

        Some(a)
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

    pub fn calc_origin_point(points: &[Point3]) -> Point3 {
        let mut temp = Point3 {x: 0.0, y: 0.0, z: 0.0};
        points.iter().for_each(|a| {
            temp.x += a.x;
            temp.y += a.y;
            temp.z += a.z;
        });
        
        Point3 { 
            x: temp.x/points.len() as f32, 
            y: temp.y/points.len() as f32, 
            z: temp.z/points.len() as f32 
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
    pub fn calc_normal(points: [Point3; 3]) -> Point3 {
        let vu = Point3 {x: points[1].x-points[0].x,y: points[1].y-points[0].y,z: points[1].z-points[0].z};
        let vv = Point3 {x: points[2].x-points[0].x,y: points[2].y-points[0].y,z: points[2].z-points[0].z};
    
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
}

impl MeshAndBuffer for Mesh {
    fn sort_triangles_in_mesh_origin_point_z_distance_from_camera(&mut self, camera: &Camera) -> Result<(), ()> {
        if self.triangles.len() == 0 {return Err(());}
        
        let check = self.triangles.iter_mut()
            .try_for_each(|a| -> Result<(), ()> { 
                if a.origin_point == None {return Err(());}
                
                a.origin_point = 
                Some(
                    Point3 {
                    x: a.origin_point.unwrap().x-camera.origin_point.x,
                    y: a.origin_point.unwrap().y-camera.origin_point.y,
                    z: a.origin_point.unwrap().z-camera.origin_point.z
                    }
                ); 
                Ok(())}
            );
        if check == Err(()) {
            return Err(());
        }

        self.triangles.sort_by(|a, b| b.origin_point.unwrap().z.partial_cmp(&a.origin_point.unwrap().z).unwrap());
        
        self.triangles.iter_mut()
            .for_each(|a| { 
                
                a.origin_point = 
                Some(
                    Point3 {
                    x: a.origin_point.unwrap().x-camera.origin_point.x,
                    y: a.origin_point.unwrap().y-camera.origin_point.y,
                    z: a.origin_point.unwrap().z-camera.origin_point.z
                    }
                ); 
            }
            );
        Ok(())
    }
}

impl Mesh {
    pub fn init(name: String, points: &[Point3], triangles: &[Triangle]) -> Mesh {
        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;

        points.iter().for_each(|a| {
            x += a.x;
            y += a.y;
            z += a.z;
        });

        let p = Point3 {x: (x/points.len() as f32), y: (y/points.len() as f32), z: (z/points.len() as f32)};
        Mesh {name: (name), points: (points.to_vec()), triangles: (triangles.to_vec()), origin_point: (p)}
    }

    pub fn set_origin_point(&mut self) {
        let mut p = Point3 { x: 0.0, y: 0.0, z: 0.0 };

        self.points.iter().for_each(|a| {
            p.x += a.x;
            p.y += a.y;
            p.z += a.z;
        });

        p.x /= self.points.len() as f32;
        p.y /= self.points.len() as f32;
        p.z /= self.points.len() as f32;
        self.origin_point = p;
    }

    pub fn rotate_x(&mut self, mut angle: f32) {
        angle *= PI/180.0;

        self.points.iter_mut().for_each(
            |a| {
            a.y -= self.origin_point.y;
            a.z -= self.origin_point.z;
            
            let y = a.y;
            let z = a.z;

            a.y = y * angle.cos() - z * angle.sin();
            a.z = y * angle.sin() + z * angle.cos();

            a.y += self.origin_point.y;
            a.z += self.origin_point.z;
        });

        self.set_normals();
    }

    pub fn rotate_y(&mut self, mut angle: f32) {
        angle *= PI/180.0;

        self.points.iter_mut().for_each(|a| {
            a.x -= self.origin_point.x;
            a.z -= self.origin_point.z;

            let x = a.x;
            let z = a.z;

            a.x = x * angle.cos() + z * angle.sin();
            a.z = -x * angle.sin() + z * angle.cos();

            a.x += self.origin_point.x;
            a.z += self.origin_point.z;
        });
        self.set_normals();
    }

    pub fn rotate_z(&mut self, mut angle: f32) {
        angle *= PI/180.0;

        self.points.iter_mut().for_each(|a| {
            a.x -= self.origin_point.x;
            a.y -= self.origin_point.y;

            let x = a.x;
            let y = a.y;

            a.x = x * angle.cos() - y * angle.sin();
            a.y = x * angle.sin() + y * angle.cos();

            a.x += self.origin_point.x;
            a.y += self.origin_point.y;
        });
        self.set_normals();
    }

    pub fn set_normals(&mut self) {
        self.triangles.iter_mut().for_each(|a|
            a.normal = Some(Triangle::calc_normal([self.points[a.point1], self.points[a.point2], self.points[a.point3]]))
        );
    }

    pub fn set_color(&mut self, color: Color) {
        self.triangles.iter_mut().for_each(|a| a.color = Some(color));
    }
}

#[derive(Debug, Clone)]
pub struct Buffer {
    amt_of_points: usize,

    pub points: Vec<Point3>,
    pub triangles: Vec<Triangle>,
}
impl Buffer {
    pub fn init() -> Buffer {
        Buffer { amt_of_points: 0, points: Vec::new(), triangles: Vec::new() }
    }

    pub fn load_mesh(&mut self, mesh: &Mesh) {
        mesh.points.iter().for_each(|a| self.points.push(*a));
        mesh.triangles.iter().for_each(|a| self.triangles.push(Triangle { point1: a.point1+self.amt_of_points, point2: a.point2+self.amt_of_points, point3: a.point3+self.amt_of_points, origin_point: a.origin_point, normal: a.normal, color: a.color }));

        self.amt_of_points += mesh.points.len();
    }

    pub fn load_meshes(&mut self, meshes: &[Mesh]) {
        meshes.iter().for_each(|a| self.load_mesh(a));
    }
}

impl<'a> MeshAndBuffer for Buffer {
    fn sort_triangles_in_mesh_origin_point_z_distance_from_camera(&mut self, camera: &Camera) -> Result<(), ()> {
        if self.triangles.len() == 0 {return Err(());}
        
        let check = self.triangles.iter_mut()
            .try_for_each(|a| -> Result<(), ()> { 
                if a.origin_point == None {return Err(());}
                
                a.origin_point = 
                Some(
                    Point3 {
                    x: a.origin_point.unwrap().x-camera.origin_point.x,
                    y: a.origin_point.unwrap().y-camera.origin_point.y,
                    z: a.origin_point.unwrap().z-camera.origin_point.z
                    }
                ); 
                Ok(())}
            );
        if check == Err(()) {
            return Err(());
        }

        self.triangles.sort_by(|a, b| b.origin_point.unwrap().z.partial_cmp(&a.origin_point.unwrap().z).unwrap());
        
        self.triangles.iter_mut()
            .for_each(|a| { 
                
                a.origin_point = 
                Some(
                    Point3 {
                    x: a.origin_point.unwrap().x-camera.origin_point.x,
                    y: a.origin_point.unwrap().y-camera.origin_point.y,
                    z: a.origin_point.unwrap().z-camera.origin_point.z
                    }
                ); 
            }
            );
        Ok(())
    }
}

pub trait Lighting {
    fn to_light(&self) -> Light;
}
impl Lighting for Point3 {
    fn to_light(&self) -> Light {
        Light::Point(
            PointLight  { 
                origin_point: Point3 { x: self.x, y: self.y, z: self.z },
                brightness: 1.0,
                color: Color { r: 255, g: 255, b: 255, a: 255 }, 
                coefficients: (1.0, 1.0, 1.0),
            }
        )
    }
}
impl Lighting for Camera {
    fn to_light(&self) -> Light {
        Light::Point(
            PointLight  { 
                origin_point: self.origin_point,
                brightness: 1.0,
                color: Color { r: 255, g: 255, b: 255, a: 255 },
                coefficients: (1.0, 1.0, 1.0), 
            }
        )
    }
}

pub trait ConvertColors {
    fn convert_255_1(color: &Color) -> (f32, f32, f32, f32) {
        (
            color.r as f32 / 255.0,
            color.g as f32 / 255.0,
            color.b as f32 / 255.0,
            color.a as f32 / 255.0,
        )
    }

    fn convert_1_255(color: &(f32, f32, f32, f32)) -> Color {
        Color { 
            r: (color.0 * 255.0) as u8,
            g: (color.1 * 255.0) as u8,
            b: (color.2 * 255.0) as u8,
            a: (color.3 * 255.0) as u8
        }
    }
}
impl ConvertColors for Color {}

pub struct PointLight {
    pub coefficients: (f32, f32, f32),

    pub origin_point: Point3,
    pub brightness: f32,
    pub color: Color,
}
impl PointLight {
    fn apply_shading(color: &Color, normal: &Point3, light: &PointLight) -> Color {
        let d = Point3::calculate_distance(&normal, &light.origin_point);

        let c = Point3::get_dot_product(normal, &light.origin_point);

        let converted_color = Color::convert_255_1(&color);
        let converted_light = Color::convert_255_1(&light.color);

        let c = ( 
            converted_color.0 * converted_light.0 * (c/d.powi(2)) * light.brightness,
            converted_color.1 * converted_light.1 * (c/d.powi(2)) * light.brightness,
            converted_color.2 * converted_light.2 * (c/d.powi(2)) * light.brightness,
            converted_color.3,
        );
        Color::convert_1_255(&c)
    }
}

pub enum Light {
    Point(PointLight),
}

impl Light {
    fn apply_shading_multiple_lights(color: &Color, normal: &Point3, lights: &[Light]) -> Color {
        let mut c = (0.0, 0.0, 0.0, 255.0);
        lights.iter().for_each(|a| {
            let temp = Color::convert_255_1(&Light::apply_shading(color, normal, a));
            c.0 += temp.0;
            c.1 += temp.1;
            c.2 += temp.2;
        });
        Color::convert_1_255(&c)
    }

    fn apply_shading(color: &Color, normal: &Point3, light: &Light) -> Color {
        match light {
            Light::Point(x) => return PointLight::apply_shading(color, normal, x)
        }
    }
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

    rotations: (f32, f32, f32),

    z_depth_buffer: Vec<f32>,

    pub forward_vec: Point3,
    pub right_vec: Point3,
    pub up_vec: Point3,
    pub look_vec: Point3,
    pub target_vec: Point3,
}

impl<'a> Camera {
    pub fn init(origin_point: Point3, width: u32, height: u32, fov: f32, near_clipping_plane: f32, far_clipping_plane: f32) -> Camera {
        Camera { 
            origin_point: (origin_point),
            width: (width),
            height: (height),
            fov: fov,
            fov_radians: (fov*PI as f32/180.0),
            near_clipping_plane: (near_clipping_plane),
            far_clipping_plane: (far_clipping_plane),
            z_depth_buffer: vec![far_clipping_plane; width as usize * height as usize],
            forward_vec: Point3::normalize_vector(Point3 { x: 0.0, y: 0.0, z: 1.0 }),
            right_vec: Point3::normalize_vector(Point3 { x: -1.0, y: 0.0, z: 0.0 }),
            up_vec: Point3::normalize_vector(Point3 { x: 0.0, y: 1.0, z: 0.0 }),
            look_vec: Point3::normalize_vector(Point3 { x: 0.0, y: 0.0, z: 1.0 }),
            target_vec: Point3::normalize_vector(Point3 { x: 0.0, y: 0.0, z: 1.0 }),
            rotations: (0.0, 0.0, 0.0)
        }
    }

    pub fn rotate_x(&mut self, angle: f32) {
        self.rotations.0 += angle;
        if self.rotations.0 == 360.0 {
            self.rotations.0 = 0.0
        }
        else if self.rotations.0 > 360.0 {
            self.rotations.0 = self.rotations.0 - 360.0
        }
    }

    pub fn rotate_y(&mut self, angle: f32) {
        self.rotations.1 += angle;
        if self.rotations.1 == 360.0 {
            self.rotations.1 = 0.0
        }
        else if self.rotations.1 > 360.0 {
            self.rotations.1 = self.rotations.1 - 360.0
        }
    }

    pub fn rotate_z(&mut self, angle: f32) {
        self.rotations.2 += angle;
        if self.rotations.2 == 360.0 {
            self.rotations.2 = 0.0
        }
        else if self.rotations.2 > 360.0 {
            self.rotations.2 = self.rotations.2 - 360.0
        }
    }
}

pub struct MeshLoader(pub Vec<Mesh>);
impl MeshLoader {
    pub fn load_obj_file(&mut self, path:&str) -> Result<(), Box<dyn Error>> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut normals: Vec<Point3> = Vec::new();
        let mut amt_of_points = 0;
        let mut counter = 0;

        for line in reader.lines() {
            let line_string = line.unwrap();
            let mut line_split = line_string.split(" ");
            let line_0 = line_split.nth(0).unwrap().to_string(); 
            let obj_len = self.0.len();

            if line_0 == "o" {
                if counter > 0 {
                    amt_of_points += self.0[obj_len-1].points.len();
                }
                self.0.push(Mesh::init(line_split.nth(0).unwrap().to_string(), &[],  &[]));
                counter += 1;
            }

            else if line_0 == "vn" {
                let x = line_split.nth(0).unwrap().to_string();
                let y = line_split.nth(0).unwrap().to_string();
                let z = line_split.nth(0).unwrap().to_string();
                normals.push(Point3 { x: x.parse()?, y: y.parse()?, z: z.parse()? });
            }
            
            else if line_0 == "v" {
                let x = line_split.nth(0).unwrap().to_string();
                let y = line_split.nth(0).unwrap().to_string();
                let z = line_split.nth(0).unwrap().to_string();
                self.0[obj_len-1].points.push(Point3 {x: x.parse()?,y: y.parse()?,z: z.parse()?});
            }

            else if line_0 == "f" {
                let mut p_1: u32;
                let mut p_2: u32;
                let mut p_3: u32;
                let mut normal: Option<Point3> = None;
                if line_string.contains("//") {
                    let part_1 = line_split.nth(0).unwrap().to_string();
                    let part_2 = line_split.nth(0).unwrap().to_string();
                    let part_3 = line_split.nth(0).unwrap().to_string();

                    let mut part_1_split = part_1.split("//");
                    let mut part_2_split = part_2.split("//");
                    let mut part_3_split = part_3.split("//");

                    p_1 = part_1_split.nth(0).unwrap().to_string().parse()?;
                    let v_1: usize = part_1_split.nth(0).unwrap().to_string().parse()?;

                    p_2 = part_2_split.nth(0).unwrap().to_string().parse()?;
                    let v_2: usize = part_2_split.nth(0).unwrap().to_string().parse()?;
                    
                    p_3 = part_3_split.nth(0).unwrap().to_string().parse()?;
                    let v_3: usize = part_3_split.nth(0).unwrap().to_string().parse()?;

                    normal = Some(normals[v_2-1]);
                }
                else {
                    p_1 = line_split.nth(0).unwrap().to_string().parse()?;
                    p_2 = line_split.nth(0).unwrap().to_string().parse()?;
                    p_3 = line_split.nth(0).unwrap().to_string().parse()?;
                }
                p_1 -= 1 + amt_of_points as u32;
                p_2 -= 1 + amt_of_points as u32;
                p_3 -= 1 + amt_of_points as u32;

                let p = Point3 {
                    x: (self.0[obj_len-1].points[p_1 as usize].x + self.0[obj_len-1].points[p_2 as usize].x + self.0[obj_len-1].points[p_2 as usize].x)/3.0,
                    y: (self.0[obj_len-1].points[p_1 as usize].y + self.0[obj_len-1].points[p_2 as usize].y + self.0[obj_len-1].points[p_2 as usize].y)/3.0,
                    z: (self.0[obj_len-1].points[p_1 as usize].z + self.0[obj_len-1].points[p_2 as usize].z + self.0[obj_len-1].points[p_2 as usize].z)/3.0
                };
                let l = Some(Color { r: 180, g: 180, b: 180, a: 255 });
                let n;
                
                if normal == None {
                    n = Triangle::calc_normal([self.0[obj_len-1].points[p_1 as usize], self.0[obj_len-1].points[p_2 as usize], self.0[obj_len-1].points[p_3 as usize]]);
                }
                else {
                    n = normal.unwrap();
                }
                self.0[obj_len-1].triangles.push(Triangle { point1: (p_1 as usize), point2: (p_2 as usize), point3: (p_3 as usize), origin_point: Some(p), color: (l), normal: Some(n) });
            }
        }  
        Ok(())
    }
}

pub trait SdlWrapper {
    fn draw_lines_w(&mut self, buffer: &mut Buffer, camera: &Camera) -> Result<(), ()>;
    
    fn draw_triangle(&mut self, points: [Point3; 3], camera: &mut Camera);
    fn draw_triangles(&mut self, buffer: &mut Buffer, camera: &mut Camera, lights: &[Light]) -> Result<(), ()>;

    fn draw_all(&mut self, buffer: &mut Buffer, camera: &Camera);

    fn check_backface_culling(first_point_: Point3, normal: Point3, camera: &Camera) -> bool {
        let mut first_point = first_point_;
        Point3::rotate_x_around_point_degrees(&mut first_point, &camera.origin_point, camera.rotations.0);
        Point3::rotate_y_around_point_degrees(&mut first_point, &camera.origin_point, camera.rotations.1);
        Point3::rotate_z_around_point_degrees(&mut first_point, &camera.origin_point, camera.rotations.2);

        normal.x * (first_point.x - camera.origin_point.x) +
        normal.y * (first_point.y - camera.origin_point.y) +
        normal.z * (first_point.z - camera.origin_point.z) > 0.0
    }
}

impl SdlWrapper for Canvas<Window> {
    fn draw_lines_w(&mut self, buffer: &mut Buffer, camera: &Camera) -> Result<(), ()> {
        let check = buffer.sort_triangles_in_mesh_origin_point_z_distance_from_camera(camera);
        if check == Err(()) {return Err(());}
        self.set_draw_color(Color::RGB(0, 0, 0));

        let check = buffer.triangles.iter().try_for_each(|a| {
            if a.normal == None {return Err(())};

            if Canvas::check_backface_culling(buffer.points[a.point1], a.normal.unwrap(), camera) {
                return Ok(());
            }
            let p1 = Point3::project_point(buffer.points[a.point1], camera);
            let p2 = Point3::project_point(buffer.points[a.point2], camera);
            let p3 = Point3::project_point(buffer.points[a.point3], camera);

            if p1 == None || p2 == None || p3 == None {
                return Ok(());
            }
            let c1 = self.draw_line(Point::new(p1.unwrap().x as i32, p1.unwrap().y as i32), Point::new(p2.unwrap().x as i32, p2.unwrap().y as i32));
            let c2 = self.draw_line(Point::new(p1.unwrap().x as i32, p1.unwrap().y as i32), Point::new(p3.unwrap().x as i32, p3.unwrap().y as i32));
            let c3 = self.draw_line(Point::new(p2.unwrap().x as i32, p2.unwrap().y as i32), Point::new(p3.unwrap().x as i32, p3.unwrap().y as i32));

            if c1 != Ok(()) || c2 != Ok(()) || c3 != Ok(()) {
                return Err(());
            }

            Ok(())
        });
        if check == Err(()) {return Err(());}
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

    fn draw_triangles(&mut self, buffer: &mut Buffer, camera: &mut Camera, lights: &[Light]) -> Result<(), ()> {
        let check = buffer.sort_triangles_in_mesh_origin_point_z_distance_from_camera(camera);
        if check == Err(()) {return Err(());}

        let check = buffer.triangles.iter().try_for_each(|a| {
            if a.normal == None {
                return Err(());
            }

            if Canvas::check_backface_culling(buffer.points[a.point1], a.normal.unwrap(), camera) {
                return Ok(());
            }

            if a.color == None {
                if !lights.is_empty() {
                    let mut color = Color {r: 180, g: 180, b: 180, a: 255};
                    color = Light::apply_shading_multiple_lights(&color, &a.normal.unwrap(), lights);
                    self.set_draw_color(color);
                }
                else {
                    //TODO: change default color
                    self.set_draw_color(Color {r: 180, g: 180, b: 180, a: 255});
                }
            }
            else {
                if !lights.is_empty() {
                    let color = Light::apply_shading_multiple_lights(&a.color.unwrap(), &a.normal.unwrap(), lights);
                    self.set_draw_color(color);
                }
                else {
                    self.set_draw_color(a.color.unwrap());
                }
            }

            let p1 = Point3::project_point(buffer.points[a.point1], camera);
            let p2 = Point3::project_point(buffer.points[a.point2], camera);
            let p3 = Point3::project_point(buffer.points[a.point3], camera);

            if p1 == None || p2 == None || p3 == None {
                return Ok(());
            }
            self.draw_triangle([p1.unwrap(), p2.unwrap(), p3.unwrap()], camera);

            Ok(())
        });
        if check == Err(()) {return Err(());}
        Ok(())
    }

    fn draw_all(&mut self, buffer: &mut Buffer, camera: &Camera) {

    }
}

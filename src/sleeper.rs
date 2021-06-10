use std::{thread, time};
pub fn sleep() {
    println!("start");
    thread::sleep(time::Duration::from_secs(2));
    println!("end");
}

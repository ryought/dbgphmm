pub mod hoge {
    #[derive(Debug)]
    struct User {
        id: String,
        active: bool,
    }
    impl User {
        fn mock() -> User {
            User {
                id: String::from("hoge"),
                active: true,
            }
        }
        fn is_active(&self) -> bool {
            self.active
        }
        fn activated(self) -> User {
            User {
                id: self.id,
                active: true,
            }
        }
        fn activate(&mut self) {
            println!("activated!");
            self.active = true;
        }
    }
    pub fn test1() {
        let mut u1 = User {
            id: String::from("hogeehoge"),
            active: false,
        };
        u1.activate();
        println!("{}", u1.is_active());
        println!("{}", u1.id);

        let u2 = User::mock();
        println!("mock {:?}", u2);
    }
}

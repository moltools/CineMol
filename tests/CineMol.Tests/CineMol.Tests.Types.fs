module CineMol.Tests.Types

open System 
open NUnit.Framework
open CineMol.Types.Fundamentals
open CineMol.Types.Geometry 
    
[<TestFixture>]
type CineMolTests () =
    
    [<Test>]
    member _.``Test two points assessed to be on opposite sides of a line`` () =
        let aPoint: Point2D = { X = -1.0; Y =  1.0 }
        let bPoint: Point2D = { X =  1.0; Y = -1.0 }
        let line = Line2D ({ X = -1.0; Y = -1.0 }, { X = 1.0; Y = 1.0 })
        Assert.False(line.SameSideOfLine aPoint bPoint)
    
    [<Test>]
    member _.``Test two points assessed to be on the same side of a line`` () =
        let aPoint: Point2D = { X = -1.0; Y = 1.0 }
        let bPoint: Point2D = { X = -2.0; Y = 2.0 }
        let line = Line2D ({ X = -1.0; Y = -1.0 }, { X = 1.0; Y = 1.0 })
        Assert.True(line.SameSideOfLine aPoint bPoint)
        
    [<Test>]
    member _.``Test calculating intersection between two intersecting lines`` () =
        let aLine = Line2D ({ X = -1.0; Y =  1.0 }, { X = 1.0; Y = -1.0 })
        let bLine = Line2D ({ X = -1.0; Y = -1.0 }, { X = 1.0; Y =  1.0 })
        match aLine.IntersectionWithLine bLine with
        | Some intersection -> Assert.AreEqual({ X = 0.0; Y = 0.0 }, intersection)
        | None -> Assert.Fail()
        
    [<Test>]
    member _.``Test calculation intersection between two non-intersecting lines with the same slope`` () =
        let aLine = Line2D ({ X = -1.0; Y = 1.0 }, { X = 1.0; Y = -1.0 })
        let bLine = Line2D ({ X = -2.0; Y = 2.0 }, { X = 2.0; Y = -2.0 })
        match aLine.IntersectionWithLine bLine with
        | Some _ -> Assert.Fail()
        | None -> Assert.Pass()
        
    [<Test>]
    member _.``Test calculation intersection between two intersecting lines with slight difference in slope`` () =
        let aLine = Line2D ({ X = -1.0; Y = 1.0 }, { X = 1.0; Y = -1.0 })
        let bLine = Line2D ({ X = -0.999; Y = 0.999 }, { X = 1.0; Y = -1.0 })
        match aLine.IntersectionWithLine bLine with
        | Some _ -> Assert.Fail()
        | None -> Assert.Pass()

    [<Test>]
    member _.``Test calculation if two intersecting circles intersect`` () =
        let aCircle = Circle2D ({ X = -1.0; Y = 0.0 }, Radius 1.5)
        let bCircle = Circle2D ({ X =  1.0; Y = 0.0 }, Radius 1.5)
        Assert.True(aCircle.IntersectsWithCircle bCircle)

    [<Test>]
    member _.``Test calculation if two touching circles intersect`` () =
        let aCircle = Circle2D ({ X = -1.0; Y = 0.0 }, Radius 1.0)
        let bCircle = Circle2D ({ X =  1.0; Y = 0.0 }, Radius 1.0)
        Assert.False(aCircle.IntersectsWithCircle bCircle)
            
    [<Test>]
    member _.``Test calculation if two non-intersecting circles intersect`` () =
        let aCircle = Circle2D ({ X = -1.0; Y = 0.0 }, Radius 0.5)
        let bCircle = Circle2D ({ X =  1.0; Y = 0.0 }, Radius 0.5)
        Assert.False(aCircle.IntersectsWithCircle bCircle)
            
    [<Test>]
    member _.``Test calculation intersection two intersecting circles`` () =
        let aCircle = Circle2D ({ X = -1.0; Y = 0.0 }, Radius 2.0)
        let bCircle = Circle2D ({ X =  1.0; Y = 0.0 }, Radius 2.0)
        let expected = { X = 0.0; Y = -Math.Sqrt(3.0) }, { X = 0.0; Y = Math.Sqrt(3.0)}
        match aCircle.IntersectionWithCircle bCircle with
        | Some result -> Assert.AreEqual(expected, result)
        | None -> Assert.Fail()
        
    [<Test>]
    member _.``Test calculation intersection two touching circles`` () =
        let aCircle = Circle2D ({ X = -1.0; Y = 0.0 }, Radius 1.0)
        let bCircle = Circle2D ({ X =  1.0; Y = 0.0 }, Radius 1.0)
        match aCircle.IntersectionWithCircle bCircle with
        | Some _ -> Assert.Fail ()
        | None -> Assert.Pass()
        
    [<Test>]
    member _.``Test calculation intersection two eclipsing circles`` () =
        let aCircle = Circle2D ({ X = 0.0; Y = 0.0 }, Radius 2.0)
        let bCircle = Circle2D ({ X = 0.0; Y = 0.0 }, Radius 1.0)
        match aCircle.IntersectionWithCircle bCircle with
        | Some _ -> Assert.Fail()
        | None -> Assert.Pass()
        
    [<Test>]
    member _.``Test calculation intersection two non-intersecting circles`` () =
        let aCircle = Circle2D ({ X = -1.0; Y = 0.0 }, Radius 0.5)
        let bCircle = Circle2D ({ X =  1.0; Y = 0.0 }, Radius 0.5)
        match aCircle.IntersectionWithCircle bCircle with
        | Some _ -> Assert.Fail()
        | None -> Assert.Pass()
        
    [<Test>]
    member _.``Test calculation if two intersecting spheres intersect`` () =
        let aSphere = Sphere ({ X = -1.0; Y = -1.0; Z = 1.0 }, Radius 2.0)
        let bSphere = Sphere ({ X =  1.0; Y =  1.0; Z = 1.0 }, Radius 2.0)
        Assert.True(aSphere.IntersectsWithSphere bSphere)
        
    [<Test>]
    member _.``Test calculation if two non-intersecting spheres intersect`` () =
        let aSphere = Sphere ({ X = -1.0; Y = -1.0; Z = -1.0 }, Radius 1.0)
        let bSphere = Sphere ({ X =  1.0; Y =  1.0; Z =  1.0 }, Radius 1.0)
        Assert.False(aSphere.IntersectsWithSphere bSphere)
            
    [<Test>]
    member _.``Test calculation intersection two intersecting spheres`` () =
        let aSphere = Sphere ({ X = -1.0; Y = -1.0; Z = 1.0 }, Radius 2.0)
        let bSphere = Sphere ({ X =  1.0; Y =  1.0; Z = 1.0 }, Radius 2.0)
        match aSphere.IntersectionWithSphere bSphere with
        | Some (Circle3D (center, radius, norm)) -> Assert.Pass() // TODO 
        | None -> Assert.Fail()
        
    [<Test>]
    member _.``Test calculation intersection two touching spheres`` () =
        let aSphere = Sphere ({ X = -1.0; Y = 0.0; Z = 1.0 }, Radius 1.0)
        let bSphere = Sphere ({ X =  1.0; Y = 0.0; Z = 1.0 }, Radius 1.0)
        match aSphere.IntersectionWithSphere bSphere with
        | Some _ -> Assert.Fail()
        | None -> Assert.Pass()
        
    [<Test>]
    member _.``Test calculation intersection two spheres with same centers when other sphere is inside this sphere`` () =
        let aSphere = Sphere ({ X = -1.0; Y = 2.0; Z = 3.0 }, Radius 2.0)
        let bSphere = Sphere ({ X = -1.0; Y = 2.0; Z = 3.0 }, Radius 1.0)
        match aSphere.IntersectionWithSphere bSphere with
        | Some _ -> Assert.Fail()
        | None -> Assert.Pass()
        
    [<Test>]
    member _.``Test calculation intersection two spheres with same centers when this sphere is inside other sphere`` () =
        let aSphere = Sphere ({ X = -1.0; Y = 2.0; Z = 3.0 }, Radius 1.0)
        let bSphere = Sphere ({ X = -1.0; Y = 2.0; Z = 3.0 }, Radius 2.0)
        match aSphere.IntersectionWithSphere bSphere with
        | Some _ -> Assert.Fail()
        | None -> Assert.Pass()
        
    [<Test>]
    member _.``Test calculation intersection two spheres with different centers when other sphere is inside this sphere`` () =
        let aSphere = Sphere ({ X = -1.0; Y = 2.0; Z = 3.0 }, Radius 2.0)
        let bSphere = Sphere ({ X = -1.0; Y = 1.5; Z = 3.0 }, Radius 0.5)
        match aSphere.IntersectionWithSphere bSphere with
        | Some _ -> Assert.Fail()
        | None -> Assert.Pass()
        
    [<Test>]
    member _.``Test calculation intersection two non-intersecting spheres`` () =
        let aSphere = Sphere ({ X = -1.0; Y = -1.0; Z = -1.0 }, Radius 1.0)
        let bSphere = Sphere ({ X =  1.0; Y =  1.0; Z =  1.0 }, Radius 1.0)
        match aSphere.IntersectionWithSphere bSphere with
        | Some _ -> Assert.Fail()
        | None -> Assert.Pass()
        
    [<Test>]
    member _.``Test calculation intersection intersecting line and sphere`` () =
        let aLine = Line3D ({ X = 0.0; Y = 0.0; Z = -10.0 }, { X = 0.0; Y = 0.0; Z = 10.0 })
        let aSphere = Sphere ({ X = 0.0; Y = 0.0; Z = 0.0 }, Radius 1.0)
        let expected = { X = 0.0; Y = 0.0; Z = -1.0 }, { X = 0.0; Y = 0.0; Z = 1.0 }
        match aLine.IntersectionWithSphere aSphere with
        | Some result -> Assert.AreEqual(expected, result)
        | None -> Assert.Fail()
        
    [<Test>]
    member _.``Test calculation intersection tangent line and sphere`` () =
        let aLine = Line3D ({ X = 1.0; Y = 0.0; Z = -1.0 }, { X = 1.0; Y = 0.0; Z = 1.0 })
        let aSphere = Sphere ({ X = 0.0; Y = 0.0; Z = 0.0 }, Radius 1.0)
        match aLine.IntersectionWithSphere aSphere with
        | Some _ -> Assert.Fail()
        | None -> Assert.Pass()
        
    [<Test>]
    member _.``Test calculation intersection non-intersecting line and sphere`` () =
        let aLine = Line3D ({ X = 1.0; Y = 0.0; Z = -1.0 }, { X = 1.0; Y = 0.0; Z = 1.0 })
        let aSphere = Sphere ({ X = -5.0; Y = -5.0; Z = -5.0 }, Radius 1.0)
        match aLine.IntersectionWithSphere aSphere with
        | Some _ -> Assert.Fail()
        | None -> Assert.Pass()
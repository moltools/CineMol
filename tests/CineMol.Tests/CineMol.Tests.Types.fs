module CineMol.Tests.Types

open NUnit.Framework
open CineMol.Types.Geometry
    
[<TestFixture>]
type CineMolTests () =
    
    [<Test>]
    member _.``Test two points assessed to be on opposite sides of a line`` () =
        let aPoint: Point2D = { X = -1.0; Y =  1.0 }
        let bPoint: Point2D = { X =  1.0; Y = -1.0 }
        let line = Line ({ X = -1.0; Y = -1.0 }, { X = 1.0; Y = 1.0 })
        Assert.False(line.SameSideOfLine aPoint bPoint)
    
    [<Test>]
    member _.``Test two points assessed to be on the same side of a line`` () =
        let aPoint: Point2D = { X = -1.0; Y = 1.0 }
        let bPoint: Point2D = { X = -2.0; Y = 2.0 }
        let line = Line ({ X = -1.0; Y = -1.0 }, { X = 1.0; Y = 1.0 })
        Assert.True(line.SameSideOfLine aPoint bPoint)
        
    [<Test>]
    member _.``Test calculating intersection between two intersecting lines`` () =
        let aLine = Line ({ X = -1.0; Y =  1.0 }, { X = 1.0; Y = -1.0 })
        let bLine = Line ({ X = -1.0; Y = -1.0 }, { X = 1.0; Y =  1.0 })
        match aLine.IntersectionWith bLine with
        | Some intersection -> Assert.AreEqual({ X = 0.0; Y = 0.0 }, intersection)
        | None -> Assert.Fail()
        
    [<Test>]
    member _.``Test calculation intersection between two non-intersecting lines with the same slope`` () =
        /// TODO
        Assert.True(false)
        
    [<Test>]
    member _.``Test calculation intersection between two intersecting lines with slight difference in slope`` () =
        /// TODO
        Assert.True(false)

    [<Test>]
    member _.``Test calculation if two intersecting circles intersect`` () =
        /// TODO
        Assert.True(false)
        
    [<Test>]
    member _.``Test calculation if two non-intersecting circles intersect`` () =
        /// TODO
        Assert.True(false)
            
    [<Test>]
    member _.``Test calculation intersection two intersecting circles`` () =
        /// TODO
        Assert.True(false)
        
    [<Test>]
    member _.``Test calculation intersection two touching circles`` () =
        /// TODO
        Assert.True(false)
        
    [<Test>]
    member _.``Test calculation intersection two eclipsing circles`` () =
        /// TODO
        Assert.True(false)
        
    [<Test>]
    member _.``Test calculation intersection two non-intersecting circles`` () =
        /// TODO
        Assert.True(false)
        
    [<Test>]
    member _.``Test calculation if two intersecting spheres intersect`` () =
        /// TODO
        Assert.True(false)
        
    [<Test>]
    member _.``Test calculation if two non-intersecting spheres intersect`` () =
        /// TODO
        Assert.True(false)
            
    [<Test>]
    member _.``Test calculation intersection two intersecting spheres`` () =
        /// TODO
        Assert.True(false)
        
    [<Test>]
    member _.``Test calculation intersection two touching spheres`` () =
        /// TODO
        Assert.True(false)
        
    [<Test>]
    member _.``Test calculation intersection two eclipsing spheres`` () =
        /// TODO
        Assert.True(false)
        
    [<Test>]
    member _.``Test calculation intersection two non-intersecting spheres`` () =
        /// TODO
        Assert.True(false)
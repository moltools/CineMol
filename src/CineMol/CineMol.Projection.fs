module CineMol.Projection

open Types.Geometry 
open Types.Svg

/// Calculate physical projection of point on canvas.
let physicalProjection camera (pov: Point3D) p =
    let pointVector = pov.FindVector p
    { X = pointVector.ProjectVector camera.Perpendicular
      Y = pointVector.ProjectVector camera.Horizon
      Z = pointVector.ProjectVector camera.Forward }

/// Apply perspective scaling factor on point on canvas.
let perspectiveProjection focalLength p =
    let scaleFactor = focalLength / p.Z
    p.Mul scaleFactor

/// Calculate projection of point from scene on canvas.
let project camera pov focalLength p =
    p |> physicalProjection camera pov |> perspectiveProjection focalLength
module CineMol.Projection

open Types.Geometry 
open Types.Svg

/// <summary>
/// Calculate physical projection of point on canvas.
/// </summary>
let physicalProjection camera (pov: Point3D) p =
    let pointVector = pov.FindVector p
    { X = pointVector.ProjectVector camera.Perpendicular
      Y = pointVector.ProjectVector camera.Horizon
      Z = pointVector.ProjectVector camera.Forward }

/// <summary>
/// Apply perspective scaling factor on point on canvas.
/// </summary>
let perspectiveProjection focalLength p =
    let scaleFactor = focalLength / p.Z
    p.Mul scaleFactor

/// <summary>
/// Calculate projection of point from scene on canvas.
/// </summary>
let project camera pov focalLength p =
    let projected = p |> physicalProjection camera pov |> perspectiveProjection focalLength
    projected.ToPoint2D()
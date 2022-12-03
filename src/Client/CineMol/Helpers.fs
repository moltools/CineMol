module Client.CineMol.Helpers

open System

/// <summary>
///     Create string representation for float.
/// </summary>
/// <param name="f">
///     Float to create string representation from.
/// </param>
/// <returns>
///     String representation for float.
/// </returns>
let floatToStr (f: float) : string = Operators.string f

/// <summary>
///     Create string representation for integer
/// </summary>
/// <param name="d">
///     Integer to create string representation from.
/// </param>
/// <returns>
///     String representation for integer.
/// </returns>
let intToStr (d: int) : string = Operators.string d

/// <summary>
///     Round float.
/// </summary>
/// <param name="n">
///     Number of digits to round float to.
/// </param>
/// <param name="f">
///     Float to round.
/// </param>
/// <returns>
///     Rounded float.
/// </returns>
let round (n: int) (f: float) = Math.Round(f, n)

/// <summary>
///     Return absolute value of float.
/// </summary>
/// <param name="f">
///     Float to calculate absolute value of.
/// </param>
/// <returns>
///     Absolute value of float.
/// </returns>
let abs f = (f ** 2.0) ** 0.5

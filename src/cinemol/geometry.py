"""
Contains classes for representing vectors and points in 2D and 3D space.
"""
import math
import random

class Vector3D:
    """
    Represents a vector in 3D space.
    """
    def __init__(self, x: float, y: float, z: float) -> None:
        """
        :param float x: x-coordinate of the vector.
        :param float y: y-coordinate of the vector.
        :param float z: z-coordinate of the vector.
        """
        self.x = x 
        self.y = y
        self.z = z

    @classmethod
    def create_random(cls) -> "Vector3D":
        """
        Creates a random vector.
        
        :return: A random vector.
        :rtype: Vector3D
        """
        return Vector3D(random.random(), random.random(), random.random())
    
    def length(self) -> float:
        """
        Calculates the length of this vector.
        
        :return: The length of this vector.
        :rtype: float
        """
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
    
    def normalize(self) -> "Vector3D":
        """
        Normalizes this vector.
        
        :return: A new vector with the same direction as this vector but with unit length.
        :rtype: Vector3D
        """
        return self.multiply(1 / self.length())
    
    def dot(self, other: "Vector3D") -> float:
        """
        Calculates the dot product of this vector with another vector.
        
        :param Vector3D other: The vector to dot with this vector.
        :return: The dot product of this vector with another vector.
        :rtype: float
        """
        return self.x * other.x + self.y * other.y + self.z * other.z
    
    def cross(self, other: "Vector3D") -> "Vector3D":
        """
        Calculates the cross product of this vector with another vector.
        
        :param Vector3D other: The vector to cross with this vector.
        :return: The cross product of this vector with another vector.
        :rtype: Vector3D
        """
        return Vector3D(
            self.y * other.z - self.z * other.y, 
            self.z * other.x - self.x * other.z, 
            self.x * other.y - self.y * other.x
        )
    
    def subtract(self, other: "Vector3D") -> "Vector3D":
        """
        Subtracts another vector from this vector.
        
        :param Vector3D other: The vector to subtract from this vector.
        :return: A new vector with the coordinates of the difference.
        :rtype: Vector3D
        """
        return Vector3D(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def multiply(self, scalar: float) -> "Vector3D":
        """
        Multiplies this vector by a scalar.
        
        :param float scalar: The scalar to multiply this vector by.
        :return: A new vector with the coordinates of the product.
        :rtype: Vector3D
        """
        return Vector3D(self.x * scalar, self.y * scalar, self.z * scalar)

class Point2D:
    """
    Represents a point in 2D space.
    """
    def __init__(self, x: float, y: float) -> None:
        """
        :param float x: x-coordinate of the point.
        :param float y: y-coordinate of the point.
        """
        self.x = x 
        self.y = y

    def subtract_point(self, other: "Point2D") -> "Point2D":
        """
        Subtracts the coordinates of another point from this point.
        
        :param Point2D other: The point to subtract from this point.
        :return: A new point with the coordinates of the difference.
        :rtype: Point2D
        """
        return Point2D(self.x - other.x, self.y - other.y)
    
    def cross(self, other: "Point2D") -> float:
        """
        Calculates the cross product of this point with another point.
        
        :param Point2D other: The point to cross with this point.
        :return: The cross product of this point with another point.
        :rtype: float
        """
        return self.x * other.y - self.y * other.x
    
class Point3D:
    """
    Represents a point in 3D space.
    """
    def __init__(self, x: float, y: float, z: float) -> None:
        """
        :param float x: x-coordinate of the point.
        :param float y: y-coordinate of the point.
        :param float z: z-coordinate of the point.
        """
        self.x = x 
        self.y = y
        self.z = z
    
    def create_vector(self, other: "Point3D") -> Vector3D:
        """
        Creates a vector from this point to another point.
        
        :param Point3D other: The other point.
        :return: A vector from this point to another point.
        :rtype: Vector3D
        """
        return Vector3D(self.x - other.x, self.y - other.y, self.z - other.z)
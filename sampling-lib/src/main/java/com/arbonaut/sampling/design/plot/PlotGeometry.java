
package com.arbonaut.sampling.design.plot;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Coordinate;


/**
 * Encapsulates creation of plot geometries.
 * 
 * Implementations may create different plot shapes such as points, circles, 
 * squares etc. 
 */
public interface PlotGeometry {
	
    
    /**
     * Returns the minimum bounding radius of the geometry.
     * Minimum bounding radius is the radius of the smallest circle, centered at
     * the plot center, that would contain all of the plot's geometry.
     * 
     * The radius may be used by sampling design inplementations to speed up 
     * some operations.                       
     */         
    public double getMBR();
    
    /**
     * Generates new sample plot geometry based on the plot centre coordinate.
     */         
    public Geometry create(Coordinate centre);
  
}


package com.arbonaut.sampling.design.cluster;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.CoordinateSequence;
import com.vividsolutions.jts.geom.impl.CoordinateArraySequence;


/**
 * I-shaped plot cluster.
 */
public class IShapedCluster implements PlotCluster {
    
    private double ydistance;
    private int numPlots;  
    
    public IShapedCluster(double ydistance, int numPlots)
    {
        this.ydistance = ydistance;
        this.numPlots = numPlots;
    }
    
    public CoordinateSequence create(Coordinate centre)
    {
        Coordinate[] coords = new Coordinate[this.numPlots];
        
        double x = centre.x;
		double y = centre.y;
        
        for (int n = 0; n != this.numPlots; n++)
        {
            coords[n] = new Coordinate(x, y);
            y+= this.ydistance;
        }
        
        return new CoordinateArraySequence(coords);
    
    }
  
}

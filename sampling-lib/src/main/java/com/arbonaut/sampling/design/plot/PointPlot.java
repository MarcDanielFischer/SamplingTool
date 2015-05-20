
package com.arbonaut.sampling.design.plot;

import org.geotools.geometry.jts.JTSFactoryFinder;

import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Coordinate;


/**
 * Constructs Point geometry for sample plots. 
 */
public class PointPlot implements PlotGeometry {
	
    private GeometryFactory gf = JTSFactoryFinder.getGeometryFactory( null );
    
    
    // Point has zero radius
    public double getMBR() { return 0; }
    
    public Geometry create(Coordinate centre)
    {
        return gf.createPoint(centre);
    }
  
}

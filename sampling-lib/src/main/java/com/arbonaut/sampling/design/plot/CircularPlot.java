
package com.arbonaut.sampling.design.plot;

import org.geotools.geometry.jts.JTSFactoryFinder;
import org.geotools.geometry.jts.GeometryBuilder;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Coordinate;


/**
 * Constructs Point geometry for sample plots. 
 */
public class CircularPlot implements PlotGeometry {
	
    private double radius;
    private GeometryBuilder builder;
    
    public CircularPlot(double radius)
    {
        this.radius = radius;
        GeometryFactory gf = JTSFactoryFinder.getGeometryFactory(null);
        this.builder = new GeometryBuilder(gf);
    }
     
    public double getMBR() { return radius; }
    
    public Geometry create(Coordinate centre)
    {
        return builder.circle(centre.x, centre.y, this.radius, 32);
    }
  
}

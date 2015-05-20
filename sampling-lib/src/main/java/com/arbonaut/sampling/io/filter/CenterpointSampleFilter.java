
package com.arbonaut.sampling.io.filter;

import java.io.IOException;
import com.arbonaut.sampling.SamplePlot;
import com.arbonaut.sampling.BadCRSException;

import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import org.opengis.referencing.NoSuchAuthorityCodeException;
import org.opengis.referencing.FactoryException;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.MathTransform;
import org.geotools.referencing.CRS;
import org.geotools.geometry.jts.JTSFactoryFinder;


/**
 * Sample writer that reduces the plot geometry to a single point placed at the
 * centre of the plot's envelope. 
 */
public class CenterpointSampleFilter implements SampleFilter {

    private GeometryFactory gf = JTSFactoryFinder.getGeometryFactory(null);

    public SamplePlot filter(SamplePlot plot) throws Exception 
    {
        Geometry geom = plot.getGeometry();
    
        // Nothing to do
        if (geom instanceof Point) return plot;

        Envelope e = geom.getEnvelopeInternal();
        if (e == null) 
            throw new Exception("Geometry has null envelope");

        Geometry newgeom = gf.createPoint(e.centre());           
        return new SamplePlot(plot, newgeom, plot.getCRS());
    }	

}

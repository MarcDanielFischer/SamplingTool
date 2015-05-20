
package com.arbonaut.sampling;

import org.opengis.referencing.crs.CoordinateReferenceSystem;
import com.vividsolutions.jts.geom.Geometry;


/**
 * This class is used to store the strata used in the sampling process as an object type of its own.
 * The reason we do it this way is that using the Stratum class we can easily associate additional properties
 * to a Geometry. 
 *  
 * The CRS property is needed so that it can be passed on to the sample Plots which can then be reprojected to 
 * LatLon using this CRS information.
 *  
 * The name property is also used to be passed on to the sample Plots and will be retrieved and written to the output file.
 */
public class Stratum {
	
	// properties
	private Geometry geometry;
	private CoordinateReferenceSystem CRS;
	private String name;


	// constructor
	public Stratum(Geometry geometry,  CoordinateReferenceSystem CRS,  String name){
		this.geometry = geometry;
		this.CRS = CRS;
		this.name = name;
	}
	
	
	// Getters and Setters
	public Geometry getGeometry() {
		return geometry;
	}



	public void setGeometry(Geometry geometry) {
		this.geometry = geometry;
	}



	public CoordinateReferenceSystem getCRS() {
		return CRS;
	}



	public void setCRS(CoordinateReferenceSystem cRS) {
		CRS = cRS;
	}



	public String getName() {
		return name;
	}


	public void setName(String name) {
		this.name = name;
	}

}

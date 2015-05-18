package production;

import java.util.ArrayList;

import org.geotools.geometry.jts.JTSFactoryFinder;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;

/**
 * This class provides static methods to construct cluster plots using cluster 
 * seed points as input data. 
 * @author daniel
 *
 */
public class Clusters {

	
	/**
	 * Grow i-shaped clusters from given seed points with a specified number of sub-plots per cluster 
	 * and a specified distance separating thesub-plots (distance in Meters).
	 * The clusters are linearly grown from south to north (the cluster seed point is the southernmost plot). 
	 * 
	 * This method calls create_I_clusters() for all input strata.
	 * 
	 * @param clusterSeedPoints points to be used as origin of each cluster (for all strata). Must be in UTM projection.
	 * @param distBetweenSubPlots distance in Meters
	 * @param numClusterSubPlots number of sub-plots per cluster
	 * @param strata this param is needed in order to check whether all generated cluster plots are inside the stratum. Must be in UTM projection. 
	 * @return
	 */
	public static ArrayList<Plot> create_I_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numClusterSubPlots, Stratum[] strata){
		ArrayList<Plot> plots = new ArrayList<Plot>();
		for (Stratum stratum : strata) {
			// filter clusterSeedPoints so that only those for each stratum remain
			ArrayList<Plot> filteredClusterSeedPoints = new ArrayList<Plot>();
			for(Plot plot : clusterSeedPoints){
				if(plot.getStratumName().equals(stratum.getName())){
					filteredClusterSeedPoints.add(plot);
				}
			}
			plots.addAll(Clusters.create_I_clusters(filteredClusterSeedPoints, distBetweenSubPlots, numClusterSubPlots, stratum));
		}
		return plots;
	}
	
	/**
	 * Grow i-shaped clusters from given seed points with a specified number of sub-plots per cluster 
	 * and a specified distance separating thesub-plots (distance in Meters).
	 * The clusters are linearly grown from south to north (the cluster seed point is the southernmost plot). 
	 * @param clusterSeedPoints points to be used as origin of each cluster. Must be in UTM projection.
	 * @param distBetweenSubPlots distance in Meters
	 * @param numClusterSubPlots number of sub-plots per cluster
	 * @param stratum this param is needed in order to check whether all generated cluster plots are inside the stratum. Must be in UTM projection. 
	 */
	public static ArrayList<Plot> create_I_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numClusterSubPlots, Stratum stratum){
		
		ArrayList<Plot> outputPlots = new ArrayList<Plot>();
		
		// iterate over seed points
		for(Plot seedPoint : clusterSeedPoints){
			
			
			// use seedPoint.plotNr as clusterNr and set new seedPoint.plotNr to 1
			int clusterNr = seedPoint.getPlotNr();
			int subPlotNr = 1;
			seedPoint.setClusterNr(clusterNr);
			seedPoint.setPlotNr(subPlotNr);
			
			
			// add each seed point to output after changing its numbering
			outputPlots.add(seedPoint);
			
			double x = seedPoint.getPoint().getCoordinate().x;
			double y = seedPoint.getPoint().getCoordinate().y;
			
			// create additional points until reaching the specified total number of plots per cluster
			for(int i = 0; i < numClusterSubPlots - 1; i++){
				
				y += distBetweenSubPlots; // build i-clusters along North-South axis (change only y coordinate)				
				
				// create Points using GeometryFactory
				GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
				Coordinate coord = new Coordinate( x, y );
				Point point = geometryFactory.createPoint( coord );

				// check if Point inside stratum, create Plot and add Plot to output ArrayList
				if(point.within(stratum.getGeometry())){
					subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
					// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
					Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
					// check if there is a weight value and add it
					double weight = seedPoint.getWeight();
					if( weight != -1.0){
						plot.setWeight(weight);
					}
					outputPlots.add(plot);
					
				}
			}
		}
		return outputPlots;
	}
	
	


	/**
	 * Grow L-shaped clusters from given seed points with a 
	 * specified number of sub-plots per cluster 
	 * and a specified distance separating the sub-plots (distance in Meters).
	 * In case of an even total number of sub-plots shaping the cluster, the vertical axis 
	 * will be one sub-plots longer than the horizontal axis (just like the real character "L"). 
	 * 
	 * This method calls create_L_clusters() for all input strata.
	 * 
	 * @param clusterSeedPoints points to be used as origin of each cluster (for all strata). Must be in UTM projection.
	 * @param distBetweenSubPlots distance in Meters
	 * @param numClusterSubPlots number of sub-plots per cluster
	 * @param orientation 1: regular "L" shape -1: upside down
	 * @param strata this param is needed in order to check whether all generated cluster plots are inside the stratum. Must be in UTM projection. 
	 * @return
	 */
	public static ArrayList<Plot> create_L_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numClusterSubPlots, int orientation,  Stratum[] strata ){
		ArrayList<Plot> plots = new ArrayList<Plot>();
		for (Stratum stratum : strata) {
			// filter clusterSeedPoints so that only those for each stratum remain
			ArrayList<Plot> filteredClusterSeedPoints = new ArrayList<Plot>();
			for(Plot plot : clusterSeedPoints){
				if(plot.getStratumName().equals(stratum.getName())){
					filteredClusterSeedPoints.add(plot);
				}
			}
			plots.addAll(Clusters.create_L_clusters(filteredClusterSeedPoints, distBetweenSubPlots, numClusterSubPlots, orientation, stratum));
		}
		return plots;
	}
	
	/**
	 * Grow L-shaped clusters from given seed points with a 
	 * specified number of sub-plots per cluster 
	 * and a specified distance separating the sub-plots (distance in Meters).
	 * In case of an even total number of sub-plots shaping the cluster, the vertical axis 
	 * will be one sub-plots longer than the horizontal axis (just like the real character "L").
	 * @param clusterSeedPoints. Must be in UTM projection.
	 * @param distBetweenSubPlots
	 * @param numClusterSubPlots
	 * @param orientation 1: regular "L" shape -1: upside down
	 * @param stratum this param is needed in order to check whether all generated cluster plots are inside the stratum. Must be in UTM projection. 
	 */
	public static ArrayList<Plot> create_L_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numClusterSubPlots, int orientation,  Stratum stratum ){
		
		// initialize output ArrayList
		ArrayList<Plot> outputPlots = new ArrayList<Plot>();

		// iterate over seed points
		for(Plot seedPoint : clusterSeedPoints){

			int clusterNr = 0;
			int subPlotNr = 0;

			if(numClusterSubPlots > 0){
				// use seedPoint.plotNr as clusterNr and set new seedPoint.plotNr to 1
				clusterNr = seedPoint.getPlotNr();
				subPlotNr = 1;
				seedPoint.setClusterNr(clusterNr);
				seedPoint.setPlotNr(subPlotNr);

				// add each seed point to output after changing its numbering
				outputPlots.add(seedPoint);

				// Determine number of plots to be created along each axis of the L
				int numPlotsVertical = 0;
				int numPlotsHorizontal = 0;

				// odd total number of Plots shaping the cluster
				if(numClusterSubPlots % 2 != 0){ 
					// in case of an odd total Plot number, both axes are of the same length
					numPlotsVertical = numPlotsHorizontal = (numClusterSubPlots-1) / 2;
				} 

				// even total number of Plots shaping the cluster
				if(numClusterSubPlots % 2 == 0){ 
					// in case of an even total Plot number, the vertical axis will be one Plot longer than the horizontal axis
					numPlotsVertical = numClusterSubPlots / 2;
					numPlotsHorizontal = (numClusterSubPlots / 2) -1;
				} 


				// extract seed point coords and use them to build the other Plots
				double x = seedPoint.getPoint().getCoordinate().x;
				double y = seedPoint.getPoint().getCoordinate().y;


				// build Plots along vertical axis
				for(int i = 0; i < numPlotsVertical; i++){

					if(orientation == 1){ // regular "L"
						y += distBetweenSubPlots; // as we build along the vertical axis, only y coords are affected
					}
					if(orientation == -1){ // upside down "L"
						y -= distBetweenSubPlots;
					}

					// create Points using GeometryFactory
					GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
					Coordinate coord = new Coordinate( x, y );
					Point point = geometryFactory.createPoint( coord );

					// check if Point inside stratum, create Plot and add Plot to output ArrayList
					if(point.within(stratum.getGeometry())){
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						subPlotNr++;
						Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						// check if there is a weight value and add it
						double weight = seedPoint.getWeight();
						if( weight != -1.0){
							plot.setWeight(weight);
						}
						outputPlots.add(plot);
					}
				}


				// reset y coord (which has been increased during vertical axis Plot construction)  so that we can start building plots along the horizontal axis
				y = seedPoint.getPoint().getCoordinate().y;

				// build plots along horizontal axis
				for(int i = 0; i < numPlotsHorizontal; i++){
					x += distBetweenSubPlots; // as we build along the vertical axis, only x coords are affected

					// create Points using GeometryFactory
					GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
					Coordinate coord = new Coordinate( x, y );
					Point point = geometryFactory.createPoint( coord );

					// check if Point inside stratum, create Plot and add Plot to output ArrayList
					if(point.within(stratum.getGeometry())){
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						subPlotNr++;
						Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						// check if there is a weight value and add it
						double weight = seedPoint.getWeight();
						if( weight != -1.0){
							plot.setWeight(weight);
						}
						outputPlots.add(plot);
					}
				}
			}
		}
		return outputPlots;
	}
	
	
	
	/**
	 * Grow H-shaped clusters from given seed points using
	 * the specified distance separating the sub-plots (distance in Meters).
	 * The horizontal line is created first, and from its end points
	 * the vertical lines are built. 
	 * 
	 * Note: plots which are located  on both horizontal AND vertical lines 
	 * are considered to belong only to vertical lines (so in case of odd number of plots in vertical  lines 
	 * the horizontal line endpoints do not belong to horizontal line but to vertical lines instead)
	 * 
	 * This method calls create_H_clusters() for all input strata.
	 * 
	 * @param clusterSeedPoints points to be used as origin of each cluster (for all strata). Must be in UTM projection.
	 * @param distBetweenSubPlots distance between sub-plots in Meters.
	 * @param numSubPlotsInVerticalLine number of sub-plots per vertical line.
	 * @param numSubPlotsinHhorizontalLine: number of sub-plots per horizontal line.
	 * @param strata this param is needed in order to check whether all generated cluster plots are inside the stratum. Must be in UTM projection. 
	 * @return
	 */
	public static ArrayList<Plot> create_H_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numSubPlotsInVerticalLine, int numSubPlotsInHorizontalLine, Stratum[] strata ){
		ArrayList<Plot> plots = new ArrayList<Plot>();
		for (Stratum stratum : strata) {
			// filter clusterSeedPoints so that only those for each stratum remain
			ArrayList<Plot> filteredClusterSeedPoints = new ArrayList<Plot>();
			for(Plot plot : clusterSeedPoints){
				if(plot.getStratumName().equals(stratum.getName())){
					filteredClusterSeedPoints.add(plot);
				}
			}
			plots.addAll(Clusters.create_H_clusters(filteredClusterSeedPoints, distBetweenSubPlots, numSubPlotsInVerticalLine, numSubPlotsInHorizontalLine, stratum));
		}
		return plots;
	}
	
	/**
	 * Grow H-shaped clusters from given seed points. 
	 * The horizontal line is created first, and from its end points
	 * the vertical lines are built. 
	 *
	 * @param clusterSeedPoints seed points from which the clusters are generated. Must be in UTM projection.
	 * @param distBetweenSubPlots distance between sub-plots in Meters.
	 * @param numSubPlotsInVerticalLine number of sub-plots per vertical line.
	 * @param numSubPlotsinHhorizontalLine: number of sub-plots per horizontal line.
	 * Note: plots which are located  on both horizontal AND vertical lines 
	 * are considered to belong only to vertical lines (so in case of odd number of plots in vertical  lines 
	 * the horizontal line endpoints do not belong to horizontal line but to vertical lines instead)
	 * @param stratum this param is needed in order to check whether all generated cluster plots are inside the stratum. Must be in UTM projection. 
	 * @return output plots
	 */
	public static ArrayList<Plot> create_H_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numSubPlotsInVerticalLine, int numSubPlotsInHorizontalLine, Stratum stratum ){
		
		// initialize output ArrayList
		ArrayList<Plot> outputPlots = new ArrayList<Plot>();

		// diese Zeile außerhalb vom for-loop, damit das Objekt nicht dauernd neu angelegt wird (--> spart Rechendaufwand???)
		GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
		
		
		
		// iterate over seed points and make a cluster out of each
		for(Plot seedPoint : clusterSeedPoints){
			
			
			// use seedPoint.plotNr as clusterNr and set index subPlotNr to 0
			int clusterNr = seedPoint.getPlotNr();
			int subPlotNr = 0;

			// extract coords from seed point in order to construct additional Plots
			double x = seedPoint.getPoint().getCoordinate().x;
			double y = seedPoint.getPoint().getCoordinate().y;
			
			// we need to record the min/max values for x so that we can later build the vertical lines using them  
			double xMax = x;
			double xMin = x;
			

			//-----------------------------------------------------------------------------------------------
			// Plots in horizontal line
			if(numSubPlotsInHorizontalLine > 0){ // only create Plots in horizontal line if there are any to create


				//-----------------------------------------------------
				// odd number of Plots in horizontal line
				//			bei ungerader Zahl von Horizontallinienpunkten:
				//			setze einen an die Stelle vom Saatpunkt und gehe dann nach links und rechts weg
				if(numSubPlotsInHorizontalLine % 2 != 0){ // if number of Plots in horizontal line is odd, the seed point itself will be the central Plot

					// Zentralpunkt
					seedPoint.setClusterNr(clusterNr);
					subPlotNr++;
					seedPoint.setPlotNr(subPlotNr);
					// add seed point to output after changing its numbering only if number of plots in horizontal line is odd
					outputPlots.add(seedPoint);


					// Punkte nach rechts (Osten) weg anlegen
					for(int i = 0; i < (numSubPlotsInHorizontalLine -1) / 2; i++){
						x += distBetweenSubPlots; // as we operate along a horizontal line here, only x values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );

						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);
						}

						xMax = x; // we need to record the min/max values for x so that we can build the vertical lines using them   
					}

					// reset x value so that we can create Plots in the other direction
					x = seedPoint.getPoint().getCoordinate().x;

					// Punkte nach links (Westen) weg anlegen
					for(int i = 0; i < (numSubPlotsInHorizontalLine -1) / 2; i++){
						x -= distBetweenSubPlots; // as we operate along a horizontal line here, only x values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );

						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);
						}

						xMin = x; // we need to record the min/max values for x so that we can build the vertical lines using them   
					}
				}


				//-----------------------------------------------------
				// even number of Plots in horizontal line
				// if number of Plots in horizontal line is even, the seed point will NOT be the central Plot
				if(numSubPlotsInHorizontalLine % 2 == 0){ 

					//-----------------------------------
					// Plots to the right of seed point

					// first Plot to the right is half the plotDist away from seedpoint
					x += distBetweenSubPlots / 2;
					Coordinate coordFirstPlotRight = new Coordinate( x, y );
					Point pointFirstPlotRight = geometryFactory.createPoint( coordFirstPlotRight );
					// check if Point inside stratum, create Plot and add Plot to output ArrayList
					if(pointFirstPlotRight.within(stratum.getGeometry())){
						// subPlotNr++; // as this is the first subPlot, we don´t want to increase subPlotNr yet 
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						Plot plot = new Plot(pointFirstPlotRight, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						// check if there is a weight value and add it
						double weight = seedPoint.getWeight();
						if( weight != -1.0){
							plot.setWeight(weight);
						}
						outputPlots.add(plot);

					}

					// restliche Punkte nach rechts (Osten) weg anlegen
					for(int i = 0; i < (numSubPlotsInHorizontalLine / 2) -1 ; i++){ // Achtung, veränderte Zählweise im For-loop
						x += distBetweenSubPlots; // as we operate along a horizontal line here, only x values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );

						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);

						}

						xMax = x; // we need to record the min/max values for x so that we can build the vertical lines using them   
					}

					// reset x value so that we can create Plots in the other direction
					x = seedPoint.getPoint().getCoordinate().x;


					//-----------------------------------
					// Plots to the left of seed point

					// first Plot to the left is half the plotDist away from seed point
					x -= distBetweenSubPlots / 2;
					Coordinate coordFirstPlotLeft = new Coordinate( x, y );
					Point pointFirstPlotLeft = geometryFactory.createPoint( coordFirstPlotLeft );
					// check if Point inside stratum, create Plot and add Plot to output ArrayList
					if(pointFirstPlotLeft.within(stratum.getGeometry())){
						subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						Plot plot = new Plot(pointFirstPlotLeft, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						// check if there is a weight value and add it
						double weight = seedPoint.getWeight();
						if( weight != -1.0){
							plot.setWeight(weight);
						}
						outputPlots.add(plot);

					}

					// Punkte nach links (Westen) weg anlegen
					for(int i = 0; i < (numSubPlotsInHorizontalLine / 2) -1 ; i++){ // Achtung, veränderte Zählweise im For-loop
						x -= distBetweenSubPlots; // as we operate along a horizontal line here, only x values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );

						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);

						}

						xMin = x; // we need to record the min/max values for x so that we can build the vertical lines using them   
					}
				}
			} else{ // if numSubPlotsinHhorizontalLine == 0
				// if there aren`t any Plots in the horizontal line,we use the seed point's x coord to be able to construct the vertical lines
				xMin = xMax = seedPoint.getPoint().getCoordinate().x; 
			}
			
			



			//-----------------------------------------------------------------------------------------------
			// Plots in vertical lines
			if(numSubPlotsInVerticalLine > 0){ // only create Plots in vertical lines if there are any to create


				// odd number of Plots in vertical lines --> central Plots are exactly in line with horizontal line Plots(it's different with even number of Plots in vertical lines)
				// Idea: create central Plots first, then create Plots to the north and after that Plots to the south
				if(numSubPlotsInVerticalLine % 2 != 0){ 
					// based on xMin and xMax from the horizontal line (line end Plots), we construct the two vertical lines 

					//---------------------------
					// vertical line on the right
					// central Plot (exactly in line with horizontal line Plots)
					x = xMax + distBetweenSubPlots; // extreme right Plot of horizontal line + distBetweenSubPlots
					Coordinate coordCentralPlotRight = new Coordinate( x, y );
					Point pointCentralPlotRight = geometryFactory.createPoint( coordCentralPlotRight );

					// check if Point inside stratum, create Plot and add Plot to output ArrayList
					if(pointCentralPlotRight.within(stratum.getGeometry())){
						subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						Plot plot = new Plot(pointCentralPlotRight, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						// check if there is a weight value and add it
						double weight = seedPoint.getWeight();
						if( weight != -1.0){
							plot.setWeight(weight);
						}
						outputPlots.add(plot);
					}

					// restliche Punkte nach oben (Norden) weg anlegen
					for(int i = 0; i < (numSubPlotsInVerticalLine -1 ) / 2 ; i++){ // Achtung, immer wieder veränderte Zählweise im For-loop
						y += distBetweenSubPlots; // as we operate along a vertical line here, only y values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );
						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);
						}
					}

					// reset y value so that we can create Plots in the other direction
					y = seedPoint.getPoint().getCoordinate().y;

					// restliche Punkte nach unten (Süden) weg anlegen
					for(int i = 0; i < (numSubPlotsInVerticalLine -1 ) / 2 ; i++){ // Achtung, immer wieder veränderte Zählweise im For-loop
						y -= distBetweenSubPlots; // as we operate along a vertical line here, only y values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );
						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);
						}
					}


					//---------------------------
					// vertical line on the left

					// reset y value so that we can create Plots afresh
					y = seedPoint.getPoint().getCoordinate().y;

					// central Plot
					x = xMin - distBetweenSubPlots; // extreme left Plot of horizontal line - distBetweenSubPlots
					Coordinate coordCentralPlotLeft = new Coordinate( x, y );
					Point pointCentralPlotLeft = geometryFactory.createPoint( coordCentralPlotLeft );
					// check if Point inside stratum, create Plot and add Plot to output ArrayList
					if(pointCentralPlotLeft.within(stratum.getGeometry())){
						subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						Plot plot = new Plot(pointCentralPlotLeft, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						// check if there is a weight value and add it
						double weight = seedPoint.getWeight();
						if( weight != -1.0){
							plot.setWeight(weight);
						}
						outputPlots.add(plot);
					}

					// restliche Punkte nach oben (Norden) weg anlegen
					for(int i = 0; i < (numSubPlotsInVerticalLine -1 ) / 2 ; i++){ // Achtung, immer wieder veränderte Zählweise im For-loop
						y += distBetweenSubPlots; // as we operate along a vertical line here, only y values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );

						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);
						}
					}

					// reset y value so that we can create Plots in the other direction
					y = seedPoint.getPoint().getCoordinate().y;

					// restliche Punkte nach unten (Süden) weg anlegen
					for(int i = 0; i < (numSubPlotsInVerticalLine -1 ) / 2 ; i++){ // Achtung, immer wieder veränderte Zählweise im For-loop
						y -= distBetweenSubPlots; // as we operate along a vertical line here, only y values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );

						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);
						}
					}


				}

				//-------------------------------------------
				// even number of Plots in vertical lines --> no Plots are exactly in line with horizontal line 
				if(numSubPlotsInVerticalLine % 2 == 0){ 
					// based on xMin and xMax from the horizontal line (line end Plots), we construct the two vertical lines 

					// the special thing here is that we have to derive the positions of the Plots next to the ones in the horizontal line 
					// using the Pythagorean theorem in order to keep up the correct distBetweenSubPlots

					//---------------------------
					// vertical line to the right
					// x value is the same for all Plots in one vertical line
					x = xMax + Math.sqrt( (distBetweenSubPlots * distBetweenSubPlots) - ( (distBetweenSubPlots /2) * (distBetweenSubPlots /2))); // here's the Pythagorean theorem

					// start Plot to the north
					// the start y coord is half the distBetweenSubPlots to the north of the horizontal line in this case 
					y = seedPoint.getPoint().getCoordinate().y + (distBetweenSubPlots / 2);
					Coordinate coordStartPlotNorthRight = new Coordinate( x, y );
					Point pointStartPlotNorthRight = geometryFactory.createPoint( coordStartPlotNorthRight );
					// check if Point inside stratum, create Plot and add Plot to output ArrayList
					if(pointStartPlotNorthRight.within(stratum.getGeometry())){
						subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						Plot plot = new Plot(pointStartPlotNorthRight, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						// check if there is a weight value and add it
						double weight = seedPoint.getWeight();
						if( weight != -1.0){
							plot.setWeight(weight);
						}
						outputPlots.add(plot);
					}

					// remaining Plots to the north
					for( int i = 0; i < (numSubPlotsInVerticalLine / 2 ) - 1; i++){ // Achtung, immer wieder veränderte Zählweise im For-loop
						y += distBetweenSubPlots; // as we operate along a vertical line here, only y values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );
						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);
						}
					}


					// start Plot to the south
					// the start y coord is half the distBetweenSubPlots to the south of the horizontal line in this case 
					y = seedPoint.getPoint().getCoordinate().y - (distBetweenSubPlots / 2);
					Coordinate coordStartPlotSouthRight = new Coordinate( x, y );
					Point pointStartPlotSouthRight = geometryFactory.createPoint( coordStartPlotSouthRight );
					// check if Point inside stratum, create Plot and add Plot to output ArrayList
					if(pointStartPlotSouthRight.within(stratum.getGeometry())){
						subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						Plot plot = new Plot(pointStartPlotSouthRight, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						// check if there is a weight value and add it
						double weight = seedPoint.getWeight();
						if( weight != -1.0){
							plot.setWeight(weight);
						}
						outputPlots.add(plot);
					}

					// remaining Plots to the south
					for( int i = 0; i < (numSubPlotsInVerticalLine / 2 ) - 1; i++){ // Achtung, immer wieder veränderte Zählweise im For-loop
						y -= distBetweenSubPlots; // as we operate along a vertical line here, only y values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );
						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);
						}
					}

					//---------------------------
					// vertical line to the left
					// x value is the same for all Plots in one vertical line
					x = xMin - Math.sqrt( (distBetweenSubPlots * distBetweenSubPlots) - ( (distBetweenSubPlots /2) * (distBetweenSubPlots /2))); // here's the Pythagorean theorem

					// start Plot to the north
					// the start y coord is half the distBetweenSubPlots to the north of the horizontal line in this case 
					y = seedPoint.getPoint().getCoordinate().y + (distBetweenSubPlots / 2);
					Coordinate coordStartPlotNorthLeft = new Coordinate( x, y );
					Point pointStartPlotNorthLeft = geometryFactory.createPoint( coordStartPlotNorthLeft );
					// check if Point inside stratum, create Plot and add Plot to output ArrayList
					if(pointStartPlotNorthLeft.within(stratum.getGeometry())){
						subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						Plot plot = new Plot(pointStartPlotNorthLeft, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						// check if there is a weight value and add it
						double weight = seedPoint.getWeight();
						if( weight != -1.0){
							plot.setWeight(weight);
						}
						outputPlots.add(plot);
					}

					// remaining Plots to the north
					for( int i = 0; i < (numSubPlotsInVerticalLine / 2 ) - 1; i++){ // Achtung, immer wieder veränderte Zählweise im For-loop
						y += distBetweenSubPlots; // as we operate along a vertical line here, only y values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );
						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);
						}
					}


					// start Plot to the south
					// the start y coord is half the distBetweenSubPlots to the south of the horizontal line in this case 
					y = seedPoint.getPoint().getCoordinate().y - (distBetweenSubPlots / 2);
					Coordinate coordStartPlotSouthLeft = new Coordinate( x, y );
					Point pointStartPlotSouthLeft = geometryFactory.createPoint( coordStartPlotSouthLeft );
					// check if Point inside stratum, create Plot and add Plot to output ArrayList
					if(pointStartPlotSouthLeft.within(stratum.getGeometry())){
						subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
						// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
						Plot plot = new Plot(pointStartPlotSouthLeft, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
						// check if there is a weight value and add it
						double weight = seedPoint.getWeight();
						if( weight != -1.0){
							plot.setWeight(weight);
						}
						outputPlots.add(plot);
					}

					// remaining Plots to the south
					for( int i = 0; i < (numSubPlotsInVerticalLine / 2 ) - 1; i++){ // Achtung, immer wieder veränderte Zählweise im For-loop
						y -= distBetweenSubPlots; // as we operate along a vertical line here, only y values are affected
						Coordinate coord = new Coordinate( x, y );
						Point point = geometryFactory.createPoint( coord );
						// check if Point inside stratum, create Plot and add Plot to output ArrayList
						if(point.within(stratum.getGeometry())){
							subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
							// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
							Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
							// check if there is a weight value and add it
							double weight = seedPoint.getWeight();
							if( weight != -1.0){
								plot.setWeight(weight);
							}
							outputPlots.add(plot);
						}
					}
				}
			}
		}
		return outputPlots;
	}
	
	
	
	/**
	 * Grow square-shaped clusters from given seed points with a specified number of sub-plots per cluster 
	 * and a specified distance separating thesub-plots (distance in Meters).
	 * 
	 * This method calls create_Square_clusters() for all input strata.
	 * 
	 * @param clusterSeedPoints points to be used as origin of each cluster (for all strata). Must be in UTM projection.
	 * @param distBetweenSubPlots distance in Meters
	 * @param numClusterSubPlots the total number of sub-plots per cluster. Note: any input number 
	 * will result in output clusters that consist of a number of sub-plots that is divisible by 4.
	 * @param strata this param is needed in order to check whether all generated cluster plots are inside the stratum. Must be in UTM projection. 
	 * @return
	 */
	public static ArrayList<Plot> create_Square_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numClusterSubPlots, Stratum[] strata){
		ArrayList<Plot> plots = new ArrayList<Plot>();
		for (Stratum stratum : strata) {
			// filter clusterSeedPoints so that only those for each stratum remain
			ArrayList<Plot> filteredClusterSeedPoints = new ArrayList<Plot>();
			for(Plot plot : clusterSeedPoints){
				if(plot.getStratumName().equals(stratum.getName())){
					filteredClusterSeedPoints.add(plot);
				}
			}
			plots.addAll(Clusters.create_Square_clusters(filteredClusterSeedPoints, distBetweenSubPlots, numClusterSubPlots, stratum));
		}
		return plots;
	}
	
	/**
	 * Grow square-shaped clusters from given seed points. 
	 * @param clusterSeedPoints seed points from which the clusters are generated. Must be in UTM projection.
	 * @param distBetweenSubPlots distance between sub-plots in Meters.
	 * @param numClusterSubPlots the total number of sub-plots per cluster. Note: any input number 
	 * will result in output clusters that consist of a number of sub-plots that is divisible by 4.
	 * @param stratum this param is needed in order to check whether all generated cluster 
	 * plots are inside the stratum. Must be in UTM projection.
	 * @return output plots
	 */
	public static ArrayList<Plot> create_Square_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, int numClusterSubPlots, Stratum stratum){
		
		// ersma nur für "runde" Werte implementieren: 4,8,16,32...
		// Formel für "runde" Werte: numClusterSubPlots % 2^(2+x) == 0 --> 2^2=4, 2^(2+1)=8, 2^(2+2)=16, ...
		// Vorgehen: Startpunkt oben links (willkürlich, oben links auf der Seite fängt man auch zu lesen an) und dann im Uhrzeigersinn durchnumerieren
		// erst wird der Startpunkt angelegt und dann die Seiten (oben, rechts, unten, links) in jeweils eigener for-Schleife
		// (Seed Point ist in der Mitte des Quadrats)

		// initialize output ArrayList
		ArrayList<Plot> outputPlots = new ArrayList<Plot>();

		// three objects to be used repeatedly inside the loop
		GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
		Coordinate coord;
		Point point;
		
		// iterate over seed points
		for(Plot seedPoint : clusterSeedPoints){

			// use seedPoint.plotNr as clusterNr
			int clusterNr = seedPoint.getPlotNr();
			int subPlotNr = 1;

			// extract coords from seedPoint in order to create cluster Plots using them
			double x = seedPoint.getPoint().getCoordinate().x;
			double y = seedPoint.getPoint().getCoordinate().y;

			// erst den Plot oben links (Staertpunkt) anlegen, danach die vier Seiten jede in ihrer eigenen Schleife
			// first Plot (oben links)
			x = x - (distBetweenSubPlots * (numClusterSubPlots / 8));
			y = y + (distBetweenSubPlots * (numClusterSubPlots / 8));

			// create Points using GeometryFactory
			coord = new Coordinate( x, y );
			point = geometryFactory.createPoint( coord );

			// check if Point inside stratum, create Plot and add Plot to output ArrayList
			if(point.within(stratum.getGeometry())){
				// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
				Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
				// check if there is a weight value and add it
				double weight = seedPoint.getWeight();
				if( weight != -1.0){
					plot.setWeight(weight);
				}
				outputPlots.add(plot);
				subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
			}


			// obere Seite: numClusterSubPlots / 4 mal die x-coord ändern und Punkte setzen
			for(int i = 0; i < (numClusterSubPlots / 4); i++){
				x = x + distBetweenSubPlots;

				// create Points using GeometryFactory
				coord = new Coordinate( x, y );
				point = geometryFactory.createPoint( coord );

				// check if Point inside stratum, create Plot and add Plot to output ArrayList
				if(point.within(stratum.getGeometry())){
					// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
					Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
					// check if there is a weight value and add it
					double weight = seedPoint.getWeight();
					if( weight != -1.0){
						plot.setWeight(weight);
					}
					outputPlots.add(plot);
					subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
				}
			}


			// rechte Seite: numClusterSubPlots / 4 mal die y-coord ändern und Punkte setzen
			for(int i = 0; i < (numClusterSubPlots / 4); i++){
				y = y - distBetweenSubPlots;

				// create Points using GeometryFactory
				coord = new Coordinate( x, y );
				point = geometryFactory.createPoint( coord );

				// check if Point inside stratum, create Plot and add Plot to output ArrayList
				if(point.within(stratum.getGeometry())){
					// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
					Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
					// check if there is a weight value and add it
					double weight = seedPoint.getWeight();
					if( weight != -1.0){
						plot.setWeight(weight);
					}
					outputPlots.add(plot);
					subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
				}
			}


			// untere Seite: numClusterSubPlots / 4 mal die x-coord ändern und Punkte setzen
			for(int i = 0; i < (numClusterSubPlots / 4); i++){
				x = x - distBetweenSubPlots;

				// create Points using GeometryFactory
				coord = new Coordinate( x, y );
				point = geometryFactory.createPoint( coord );

				// check if Point inside stratum, create Plot and add Plot to output ArrayList
				if(point.within(stratum.getGeometry())){
					// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
					Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
					// check if there is a weight value and add it
					double weight = seedPoint.getWeight();
					if( weight != -1.0){
						plot.setWeight(weight);
					}
					outputPlots.add(plot);
					subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
				}
			}


			// linke Seite: (numClusterSubPlots / 4) -1 mal die y-coord ändern und Punkte setzen (einmal weniger, weil der Startpunkt oben links ja schon da ist)
			for(int i = 0; i < (numClusterSubPlots / 4) -1; i++){
				y = y + distBetweenSubPlots;

				// create Points using GeometryFactory
				coord = new Coordinate( x, y );
				point = geometryFactory.createPoint( coord );

				// check if Point inside stratum, create Plot and add Plot to output ArrayList
				if(point.within(stratum.getGeometry())){
					// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
					Plot plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
					// check if there is a weight value and add it
					double weight = seedPoint.getWeight();
					if( weight != -1.0){
						plot.setWeight(weight);
					}
					outputPlots.add(plot);
					subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
				}
			}

		}
		return outputPlots;
	}
	
	
	
	/**
	 * Grow rotated square clusters from given seed points using 
	 * a specified distance separating the sub-plots (distance in Meters).
	 * As for now, this method only creates clusters containing 4 sub-Plots.
	 * 
	 * This method calls create_rotated_Square_clusters() for all input strata.
	 * 
	 * @param clusterSeedPoints points to be used as origin of each cluster (for all strata). Must be in UTM projection.
	 * @param distBetweenSubPlots distance in Meters
	 * @param strata this param is needed in order to check whether all generated cluster plots are inside the stratum. Must be in UTM projection. 
	 * @return
	 */
	public static ArrayList<Plot> create_rotated_Square_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, Stratum[] strata){
		ArrayList<Plot> plots = new ArrayList<Plot>();
		for (Stratum stratum : strata) {
			// filter clusterSeedPoints so that only those for each stratum remain
			ArrayList<Plot> filteredClusterSeedPoints = new ArrayList<Plot>();
			for(Plot plot : clusterSeedPoints){
				if(plot.getStratumName().equals(stratum.getName())){
					filteredClusterSeedPoints.add(plot);
				}
			}
			plots.addAll(Clusters.create_rotated_Square_clusters(filteredClusterSeedPoints, distBetweenSubPlots, stratum));
		}
		return plots;
	}
	
	/**
	 * Grow rotated square clusters from given seed points. As for now, this method only creates clusters containing 4 sub-Plots.
	 * @param clusterSeedPoints seed points from which the clusters are generated. Must be in UTM projection.
	 * @param distBetweenSubPlots distance between sub-plots in Meters.
	 * @param stratum this param is needed in order to check whether all generated cluster 
	 * plots are inside the stratum. Must be in UTM projection.
	 * @return  output plots
	 */
	public static ArrayList<Plot> create_rotated_Square_clusters(ArrayList<Plot> clusterSeedPoints, int distBetweenSubPlots, Stratum stratum){

			// initialize output ArrayList
			ArrayList<Plot> outputPlots = new ArrayList<Plot>();

			// three objects to be used repeatedly inside the loop
			GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory( null );
			Coordinate coord;
			Point point;
			Plot plot;
			
			// calculate Distance between cluster Plots and seed point (applying Pythagorean theorem)
			double distFromSeedPoint = distBetweenSubPlots / Math.sqrt(2.0);
			
			// iterate over seed points
			for(Plot seedPoint : clusterSeedPoints){

				// use seedPoint.plotNr as clusterNr
				int clusterNr = seedPoint.getPlotNr();
				int subPlotNr = 1;

				// extract coords from seedPoint in order to create cluster Plots using them
				double x = seedPoint.getPoint().getCoordinate().x;
				double y = seedPoint.getPoint().getCoordinate().y;

				// create cluster Plots clockwise starting from the northern Plot 
				
				// 1. Plot
				// create Point using GeometryFactory
				coord = new Coordinate( x, (y + distFromSeedPoint) );// only Y coord affected as we move north
				point = geometryFactory.createPoint( coord );

				// check if Point inside stratum, create Plot and add Plot to output ArrayList
				if(point.within(stratum.getGeometry())){
					// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
					plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
					// check if there is a weight value and add it
					double weight = seedPoint.getWeight();
					if( weight != -1.0){
						plot.setWeight(weight);
					}
					outputPlots.add(plot);
					subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
				}
				
				// 2. Plot
				// create Point using GeometryFactory
				coord = new Coordinate( (x + distFromSeedPoint ), y );// only X coord affected as we move east
				point = geometryFactory.createPoint( coord );

				// check if Point inside stratum, create Plot and add Plot to output ArrayList
				if(point.within(stratum.getGeometry())){
					// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
					plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
					// check if there is a weight value and add it
					double weight = seedPoint.getWeight();
					if( weight != -1.0){
						plot.setWeight(weight);
					}
					outputPlots.add(plot);
					subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
				}

				// 3. Plot
				// create Point using GeometryFactory
				coord = new Coordinate( x, (y - distFromSeedPoint) );// only Y coord affected as we move south
				point = geometryFactory.createPoint( coord );

				// check if Point inside stratum, create Plot and add Plot to output ArrayList
				if(point.within(stratum.getGeometry())){
					// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
					plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
					// check if there is a weight value and add it
					double weight = seedPoint.getWeight();
					if( weight != -1.0){
						plot.setWeight(weight);
					}
					outputPlots.add(plot);
					subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
				}

				// 4. Plot
				// create Point using GeometryFactory
				coord = new Coordinate( (x - distFromSeedPoint ), y );// only Y coord affected as we move west
				point = geometryFactory.createPoint( coord );

				// check if Point inside stratum, create Plot and add Plot to output ArrayList
				if(point.within(stratum.getGeometry())){
					// a plot contains -aside from the Point object as a property - the name of the stratum it is located in and CRS information
					plot = new Plot(point, stratum.getName(), stratum.getCRS(), subPlotNr, clusterNr);
					// check if there is a weight value and add it
					double weight = seedPoint.getWeight();
					if( weight != -1.0){
						plot.setWeight(weight);
					}
					outputPlots.add(plot);
					subPlotNr++; // increase subPlotNr only if generated point falls within Geometry and is succesfully added to output 
				}

			}
			return outputPlots;
	}
}

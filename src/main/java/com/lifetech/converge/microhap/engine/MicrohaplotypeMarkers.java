package com.lifetech.converge.microhap.engine;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;

import com.lifetech.converge.microhap.domain.AeComparator;
import com.lifetech.converge.microhap.domain.Marker;
import com.lifetech.converge.microhap.domain.OrderBy;
import com.lifetech.converge.microhap.domain.SortBy;

/**
 * microhaplotype marker information
 */


/**
 * Container for Marker instances
 */
public class MicrohaplotypeMarkers implements Iterable<Marker>{
    //list of instances of Marker
	private List<Marker> markers = new ArrayList<Marker>();
	
	public MicrohaplotypeMarkers(final File inputFile) throws FileNotFoundException{
		Scanner in = new Scanner(inputFile);
		// read line by line
		while(in.hasNextLine()){
			String line = in.nextLine();
			if( !(line.startsWith("#") || line.startsWith("track") || line.startsWith("locus_nickname")) ){
				String[] cols = line.split("\t");
				Marker marker = new Marker.Builder().
						markerName(cols[0]).markerId(cols[1]).
						probOfDetMix(Float.valueOf(cols[2])).avgAe(Float.valueOf(cols[3])).
						detMixRank(Integer.valueOf(cols[4])).AeRank(Integer.valueOf(cols[5])).
			            chr(cols[6]).pMostPos(Integer.valueOf(cols[7])).snps(cols[8]).build();
				markers.add(marker);
			}
		}
		in.close();		
	}
	
	
	private class MicrohaplotypeIterator implements Iterator<Marker> {
	
		private int index = 0;

		public boolean hasNext() {
			return index < markers.size();
		}

		
		public Marker next() {
			return markers.get(index++);
		}


		public void remove() {
			// TODO Auto-generated method stub
			
		}
	}
	
	public void setSortBy(SortBy sortBy, OrderBy orderBy){
		
		switch(sortBy){
		case Ae:
			switch(orderBy){
			case ASC:
				Collections.sort(markers,AeComparator.ASC);
				break;
			case DESC:
				Collections.sort(markers,AeComparator.DESC);
				break;
			default:
				break;
			}
			break;
		//TODO: implement ProbOfDetMix if it is necessary
		case ProbOfDetMix:
			break;
		default:
			break;
		}
	}
	
	public Iterator<Marker> iterator() {
		
		return new MicrohaplotypeIterator();
	}
	
	public int size(){
		return markers.size();
	}
}

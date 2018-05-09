package com.lifetech.converge.microhap.engine;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.lifetech.converge.microhap.domain.HotSpot;

public class MicrohaplotypeHotSpots {
	
	//key: markerId and value: a set of instances of Hotspot
	private Map<String,SortedSet<HotSpot>> markers = new HashMap<String,SortedSet<HotSpot>>();

	
	public MicrohaplotypeHotSpots(File bedFile) throws FileNotFoundException{

		Scanner in = new Scanner(bedFile);
		while(in.hasNextLine()){
			String line = in.nextLine();
			if( !(line.startsWith("#") || line.startsWith("track") || line.startsWith(" ")) ){
				String [] cols = line.split("\t");
				String rsId = cols[3].split("_")[1], microhapId = cols[7], chr = cols[0];
				int pos = Integer.valueOf(cols[2]);
				
				if( markers.containsKey(microhapId) ){					       //microhap markerID exists
					SortedSet<HotSpot> hotspots = markers.get(microhapId);
					hotspots.add(new HotSpot(rsId,chr,pos));
					
				}else{
					HotSpot hs = new HotSpot(rsId,chr,pos);
					SortedSet<HotSpot> rs = new TreeSet<HotSpot>();
					rs.add(hs);
					markers.put(microhapId, rs);
				}
			}			
		}
		in.close();
		
	}
	
	
	public SortedSet<HotSpot> get(String key){

		return markers.get(key);
		
	}
	
	public Set<String> keySet(){
		return markers.keySet();
	}
 	
	public Set<Entry<String, SortedSet<HotSpot>>> entrySet(){
		return markers.entrySet();
	}
	
	
}

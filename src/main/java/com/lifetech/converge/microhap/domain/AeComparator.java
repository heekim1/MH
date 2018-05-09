package com.lifetech.converge.microhap.domain;

import java.util.Comparator;

public class AeComparator implements Comparator<Marker>{
	public static final AeComparator ASC = new AeComparator(1);
	public static final AeComparator DESC = new AeComparator(-1);
    
    private int order;
    
    private AeComparator(int order){
    	this.order = order;
    }
    

	public int compare(Marker o1, Marker o2) {
		
		if(o1.getAvgAe() < o2.getAvgAe()){
			return -1 * order;
		}else if(o1.getAvgAe() == o2.getAvgAe()){
			return 0;
		}else{
			return 1 * order;
		}
		
	}

}

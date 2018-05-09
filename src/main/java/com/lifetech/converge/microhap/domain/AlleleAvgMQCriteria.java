package com.lifetech.converge.microhap.domain;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class AlleleAvgMQCriteria implements Criteria<Allele> {
	private int MIN;
	private int MAX;
	
	public AlleleAvgMQCriteria (){
		this.MIN = 0;
		this.MAX = 254;

	}
	
	public List<Allele> meetCriteria(List<Allele> alleles) {
		List<Allele> alleleList = new ArrayList<Allele>();
		for(Allele allele : alleles){
			if(allele.getAvgMQ() >= MIN && allele.getAvgMQ() <= MAX){
				alleleList.add(allele);
			}
		}
		
		return alleleList;
	}
	
	public Map<String, Allele> meetCriteria(Map<String, Allele> alleles) {
		Map<String,Allele> alleleMap = new HashMap<String,Allele>();
		
		for(String mhSeq: alleles.keySet()){
			Allele allele = alleles.get(mhSeq);

			if(allele.getAvgMQ() >= MIN && allele.getAvgMQ() <= MAX){
				alleleMap.put(mhSeq, allele);
			}
		}
		return alleleMap;
	}
	
	public int getMIN() {
		return this.MIN;
	}

	public void setMIN(int MIN) {
		this.MIN = MIN;
	}

	public int getMAX() {
		return this.MAX;
	}

	public void setMAX(int MAX) {
		this.MAX = MAX;
	}

	@Override
	public String toString() {
		return "AlleleAvgMQCriteria [MIN=" + MIN + ", MAX=" + MAX + "]";
	}
	
	
}

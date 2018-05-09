package com.lifetech.converge.microhap.domain;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class OrCriteria implements Criteria<Allele> {

	private Criteria<Allele> criteria;
	private List<Criteria<Allele>> criteriaList;
	
	public OrCriteria(Criteria<Allele> criteria, Criteria<Allele> otherCriteria){
        this.criteria = criteria;
		this.criteriaList = new ArrayList<Criteria<Allele>>();
		this.criteriaList.add(otherCriteria);

	}
	
	public Criteria<Allele> Or(Criteria<Allele> otherCriteria){
		this.criteriaList.add(otherCriteria);
		return this;
	}
	
	public List<Allele> meetCriteria(List<Allele> items) {
		
		List<Allele> alleleList = criteria.meetCriteria(items);
		for(Criteria<Allele> criteriaItem: criteriaList){
			List<Allele> otherAlleleList = criteriaItem.meetCriteria(alleleList);
			for(Allele allele: otherAlleleList){
				if(!alleleList.contains(allele)){
					alleleList.add(allele);
				}
			}
		}
		return alleleList;
	}

	public Map<String, Allele> meetCriteria(Map<String, Allele> entries) {
		Map<String,Allele> alleleMap = criteria.meetCriteria(entries);
		for(Criteria<Allele> otherCriteria: criteriaList){
			Map<String,Allele> otherAlleleMap = otherCriteria.meetCriteria(alleleMap);
			for(String mhSeq: otherAlleleMap.keySet()){
				if(!alleleMap.containsKey(mhSeq)){
					alleleMap.put(mhSeq, otherAlleleMap.get(mhSeq));
				}
			}
		}
		return alleleMap;
	}

	@Override
	public String toString() {
		return "OrCriteria [criteria=" + criteria + ", criteriaList=" + criteriaList + "]";
	}
}

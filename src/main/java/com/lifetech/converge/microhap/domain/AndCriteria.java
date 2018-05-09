package com.lifetech.converge.microhap.domain;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class AndCriteria implements Criteria<Allele> {

	private Criteria<Allele> criteria;
	private List<Criteria<Allele>> criteriaList;
	
	public AndCriteria(Criteria<Allele> criteria, Criteria<Allele> otherCriteria){
        this.criteria = criteria;
		this.criteriaList = new ArrayList<Criteria<Allele>>();
		this.criteriaList.add(otherCriteria);

	}
	
	public AndCriteria And(Criteria<Allele> otherCriteria){
		this.criteriaList.add(otherCriteria);
		return this;
	}
	
	public List<Allele> meetCriteria(List<Allele> items) {
		
		List<Allele> alleleList = criteria.meetCriteria(items);
		for(Criteria<Allele> criteriaItem: criteriaList){
			alleleList = criteriaItem.meetCriteria(alleleList);
		}
		return alleleList;
	}

	public Map<String, Allele> meetCriteria(Map<String, Allele> entries) {
		
		Map<String,Allele> alleleMap = criteria.meetCriteria(entries);
		for(Criteria<Allele> otherCriteria: criteriaList){
			alleleMap = otherCriteria.meetCriteria(alleleMap);
		}
		return alleleMap;
	}

	@Override
	public String toString() {
		return "AndCriteria [criteria=" + criteria + ", criteriaList=" + criteriaList + "]";
	}

}

package com.lifetech.converge.microhap.domain;

import java.util.List;
import java.util.Map;

public interface Criteria<T> {
	public List<T> meetCriteria(List<T> items);
	public Map<String,Allele> meetCriteria(Map<String,Allele> entries);
}

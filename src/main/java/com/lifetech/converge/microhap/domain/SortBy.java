package com.lifetech.converge.microhap.domain;

public enum SortBy {
	Ae("Ae"), ProbOfDetMix("ProbOfDetMix");
	private String name;
	
	SortBy(String name){
		this.name = name;
	}
	
	public String getName() {
        return name;
    }
}

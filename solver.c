/******* solver.c *******//*
Copyright (C) 2014 Ivan Markin
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>. */

#include <stdio.h>
#include <math.h>
#include <libnumeric.h>
#include <libconfig.h>

int main(int argc, char *argv[]){

	Complex * psi; // Wave function snapshot
	mesh space = mesh_default;
	mesh time = mesh_default;
	double p[5] = {0., 0., 0., 0., 0.}; // Packet initial parameters
	
	/// Config loading 
	config_t config;
	config_init(&config);
	if (config_read_file(&config, (argc==1) ? "libnumeric.cfg" : argv[1] )==CONFIG_FALSE) {
		printf("Configuration problem: %s\n",config_error_text(&config));		
		return -1;
	}
	
	if (	
	config_lookup_float(&config,"mesh.space.left",&space.left) *
	config_lookup_float(&config,"mesh.space.right",&space.right) *	
	config_lookup_int(&config,"mesh.space.points",&space.points) *	
	
	config_lookup_float(&config,"mesh.time.left",&time.left) *	
	config_lookup_float(&config,"mesh.time.right",&time.right) *	
	config_lookup_int(&config,"mesh.time.points",&time.points) *
	
	config_lookup_float(&config,"init.gauss_packet.mean",&p[2]) *	
	config_lookup_float(&config,"init.gauss_packet.sigma",&p[3]) *	
	config_lookup_float(&config,"init.gauss_packet.p0",&p[4]) 
	== 0
	) {
		printf("Config parameters problem! Try heck points in floats!\n");
		return -1;
	}
	config_destroy(&config);
	// Config loaded


	recalc(&space);
	recalc(&time);


	psi = (Complex *) malloc(space.points*sizeof(Complex)); //Include into recalc? 
	mesh_flush(psi,&space,&packet_framework, 0, p);
	/* Choose solving method */

	//solve_shrodinger_sweep(&space, &time, psi);
	//solve_poisson_fft(&space, psi, psi);	
	solve_ICN(&space, &time, psi);
	
	/* * * * */

	printf("Done, writing result...\n");
	write_power_to_csv("results/psi.csv",psi,&space);
	
	// Rename to destroy!
	destruct(&space);
	destruct(&time);
	return 0;
}


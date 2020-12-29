#include <stdio.h>
#include <stdlib.h>
#include "init.h"
#include "read_input.h"
#include "io_config.h"

void write_restart(char *restart_file){
  
  FILE *fid = NULL;
  fid = fopen(restart_file, "wb");
  if (fid == NULL) {
    perror("Error while creating restart file\n");
    exit(EXIT_FAILURE);
  }

  fwrite(&in, sizeof(struct input), 1, fid);
  fwrite(&sim_box_info, sizeof(struct box_info), 1, fid);
  fwrite(&part_info, sizeof(struct p_info), 1, fid);
  fwrite(part, sizeof(double), part_info.NN*4, fid);
  fclose(fid);

}


void read_restart(char *restart_file){
  
  FILE *fid = NULL;
  fid = fopen(restart_file, "rb");
  if (fid == NULL) {
    perror("Error while creating restart file\n");
    exit(EXIT_FAILURE);
  }

  struct input in_read;
  struct box_info sim_box_info_read;
  struct p_info part_info_read;
  double (*part_read)[4]; 
  part_read = malloc(part_info.NN * sizeof(*part_read));
  if (part_read == NULL){
    printf("ERROR: Failed particle allocation\n");
    exit(EXIT_FAILURE);
  }
  
  fread(&in_read, sizeof(struct input), 1, fid);
  fread(&sim_box_info_read, sizeof(struct box_info), 1, fid);
  fread(&part_info_read, sizeof(struct p_info), 1, fid);
  fread(part_read, sizeof(double), part_info.NN*4, fid);
  fclose(fid);
  
  printf("-------------------------------------\n");
  printf("%.8e %.8e\n",in.rho,in_read.rho);
  printf("%d %d\n",in.nx,in_read.nx);
  printf("%.8e %.8e\n",in.press,in_read.press);
  printf("%lu %lu\n",in.seed,in_read.seed);
  printf("%d %d\n",in.mu_insertions,in_read.mu_insertions);
  printf("%.8e %.8e\n",in.cavity_pcav,in_read.cavity_pcav);
  printf("%.8e %.8e\n",in.cavity_out_dr,in_read.cavity_out_dr);
  printf("%d %d\n",in.cavity_sample_int,in_read.cavity_sample_int);
  printf("-------------------------------------\n");
  for (int ii=0; ii<part_info.NN; ii++){
    //printf("%.8e %.8e\n",part_read[ii][1],part[ii][1]);
    if(part_read[ii][0]!=part[ii][0] || part_read[ii][1]!=part[ii][1] || part_read[ii][2]!=part[ii][2] || part_read[ii][2]!=part[ii][2]){
      printf("!!!\n");
    }
  }
  
  free(part_read);
}

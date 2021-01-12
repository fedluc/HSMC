#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <argp.h>
#include "read_input.h"
#include "nvt.h"
#include "npt.h"
//#include "cavity_nvt.h"

// ----------------------------------------
// ------ Set-up command line parser ------
// ----------------------------------------

// Documentation
static char doc[] =
  "hsmc performs Monte Carlo simulations of "
  "simple hard-sphere systems. The program is "
  "run as ./hsmc -i IN_FILE where IN_FILE is a "
  "file describing the parameters of the simulation";

// Optional arguments
static struct argp_option options[] = {
  {"input",    'i', "IN_FILE", 0,
   "Input read from IN_FILE instead of from  in.dat"}, 
  {"example",  'e', 0, 0, "Print example of input file on screen"},
  {"output",   'o', "OUT_FILE", 0,
   "Output to OUT_FILE instead of standard output" },
  { 0 }
};

// Structure to communicate between main and parser
struct arguments
{
  char *input_file;
  char *output_file;
  bool example;
};


// Single option parser
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{

  struct arguments *arguments = state->input;

  switch (key)
    {

    case 'i':
      arguments->input_file = arg;
      break;
    case 'o':
      arguments->output_file = arg;
      break;
    case 'e':
      arguments->example = true;
      break;

    case ARGP_KEY_ARG:
      if (state->arg_num > 0) // Too many necessary arguments
        argp_usage (state);

      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static struct argp argp = { options, parse_opt, 0, doc };

// ----------------------------------------
// ---------------- Main ------------------
// ----------------------------------------


int main (int argc, char **argv){

  struct arguments arguments;
  bool run=true;

  // Default values for optional arguments
  arguments.input_file  = "in.dat";
  arguments.output_file = "";
  arguments.example = false;

  // Parse command line
  argp_parse (&argp, argc, argv, 0, 0, &arguments);

  // Print example on screen
  if (arguments.example){
    print_example();
    run = false;
  }

  // Pipe output to file
  if (strcmp(arguments.output_file,"") != 0){
    if (!freopen(arguments.output_file, "w", stdout))
    {
      printf("Failed to pipe output to file %s\n",
  	     arguments.output_file);
      exit(EXIT_FAILURE);
    }
  }

  
  // Run simulation
  if (run){

    // Read input
    read_input_file(arguments.input_file);

    // Select which type of simulation to run
    if (in.press > 0){
      // NpT simulation
      hs_npt();
    }
    else if (in.cavity_pcav > 0){
      // NVT cavity simulation
      //cavity_hs_nvt();
    }
    else {
      // NVT simulation
      hs_nvt();
    }

    printf("Simulation complete!\n");

  }

  return 0;

}


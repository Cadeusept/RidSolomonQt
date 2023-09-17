#pragma once

#ifndef RID_SOLOMON_H
#define RID_SOLOMON_H

#include <iostream>
// #include <fstream>

int code_file(FILE *input_file, const char *input_filename, FILE *archived_file);

int decode_file(FILE *archived_file, FILE *output_file);

#endif

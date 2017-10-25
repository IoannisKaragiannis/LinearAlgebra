/*============================================================================
 * Name         : mylog.h is a logger for storing exceptions
 *                thrown from myLinearAlgebra library.
 * Version      : 1.0.0, 11 Sep 2017
 *
 * Copyright (c) 2017 Ioannis Karagiannis
 * All rights reserved

 * This file is part of the LinearAlgebra library.

 * LinearAlgebra is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * You are free to use this library under the terms of the GNU General
 * Public License, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with LinearAlgebra.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact info: https://www.linkedin.com/in/ioannis-karagiannis-7129394a/
 * 				ioanniskaragiannis1987@gmail.com
================================================================================*/

#ifndef MYLOG_H_
#define MYLOG_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdexcept>    // for exception, runtime_error, out_of_range
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>

// This logger will only work for Linux OS
// If you're working on Windows you should
// use the CreateDirectory() API and modify
// the following functions  accordingly
#ifdef __linux__
#define LOG_FOLDER "/tmp/LinearAlgebra"
#define LOG_ERROR_FILE "/tmp/LinearAlgebra/error.txt"
#define LOG_FILE "/tmp/LinearAlgebra/log.txt"
#define WARNING_FILE "/tmp/LinearAlgebra/warning.txt"
#endif


inline void create_directory(const std::string& folder_name){
	struct stat st = {0};
	if (stat(folder_name.c_str(), &st) == -1) {
		mkdir(folder_name.c_str(), 0700);
	}
}

inline bool file_exists (const std::string& name) {
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	} else {
		return false;
	}
}

inline void clear_file(const std::string& file){
	if(file_exists(file.c_str())){
		if( remove( file.c_str() ) != 0 ){
			perror( "Error deleting error-log file file \n");
		}
	}
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
inline const std::string my_currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
	return buf;
}

inline void log_error(const char* message){

	// Create directory if it doesn't exist.
	static uint16_t c = 0;
	if(c == 0){
		create_directory(LOG_FOLDER);
		c++;
	}

	// Concatenate time with error message
	std::string complete_message = std::string("[") + my_currentDateTime() + std::string("]") + std::string(message) + std::string("\n");
	message = complete_message.c_str();

	FILE * pFile;
	pFile = fopen(LOG_ERROR_FILE, "a+");
	if (pFile != NULL){
		fputs (message, pFile);
		fclose (pFile);
	}else{
		std::cerr << "Failed to open error file." << std::endl;
	}
}

inline void log(const char* message){

	// Create directory if it doesn't exist.
	static uint16_t c = 0;
	if(c == 0){
		create_directory(LOG_FOLDER);
		c++;
	}

	// Concatenate time with error message
	std::string complete_message = std::string("[") + my_currentDateTime() + std::string("]") + std::string(message) + std::string("\n");
	message = complete_message.c_str();

	FILE * pFile;
	pFile = fopen(LOG_FILE, "a+");
	if (pFile != NULL){
		fputs (message, pFile);
		fclose (pFile);
	}else{
		std::cerr << "Failed to open log file." << std::endl;
	}
}

inline void warning(const char* message){

	// Create directory if it doesn't exist.
	static uint16_t c = 0;
	if(c == 0){
		create_directory(LOG_FOLDER);
		c++;
	}

	// Concatenate time with error message
	std::string complete_message = std::string("[") + my_currentDateTime() + std::string("]") + std::string(message) + std::string("\n");
	message = complete_message.c_str();

	FILE * pFile;
	pFile = fopen(WARNING_FILE, "a+");
	if (pFile != NULL){
		fputs (message, pFile);
		fclose (pFile);
	}else{
		std::cerr << "Failed to open warning file." << std::endl;
	}
}


#endif /* MYLOG_H_ */

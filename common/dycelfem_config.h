// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************************************************
 *  dycelfem_config.h
 *
 *  Minimal key = value configuration reader for DyCelFEM.
 *
 *  A run is described by a single plain-text file, so that nothing about a
 *  simulation (model files, length, view, output) requires recompiling or
 *  hand-editing the C++ any more.
 *
 *      # comment
 *      input = init7140.txt...
 *
 *      [sbml]                  # section prefixes the key: sbml.common
 *      common = mcommon.xml
 *
 *      [run]
 *      steps = 300
 *
 *  Lookup order for every setting is
 *
 *      1. environment variable   DYCELFEM_<KEY>   ('.' -> '_', upper-cased)
 *      2. the config file
 *      3. the built-in default
 *
 *  so the environment always wins and can override a config file for one-off
 *  runs, and a program with no config file at all behaves exactly as before.
 ****************************************************************************/

#ifndef DYCELFEM_CONFIG_H
#define DYCELFEM_CONFIG_H

#include <string>

namespace dycfg
{
	/* Parse path. Returns false if it cannot be read. */
	bool load(const std::string &path);

	/* True once load() has succeeded. */
	bool loaded();

	/* Path of the loaded config, or "" . */
	std::string path();

	/* Heuristic: does this file look like a config file rather than a cell
	   database? True when it parses as key=value and defines "input". Lets
	   the program accept either a config or a legacy database on argv[1]. */
	bool looksLikeConfig(const std::string &path);

	/* True if the key is present in the file or the environment. */
	bool has(const std::string &key);

	std::string getStr (const std::string &key, const std::string &def);
	int         getInt (const std::string &key, int def);
	double      getNum (const std::string &key, double def);
	/* Accepts 1/0, true/false, yes/no, on/off (case-insensitive). */
	bool        getBool(const std::string &key, bool def);

	/* Resolve a filename relative to the directory holding the config file,
	   so a config can be run from anywhere. Absolute paths pass through. */
	std::string resolve(const std::string &filename);

	/* Print every key that was actually read, with its source. */
	void dumpEffective();
}

#endif

// SPDX-License-Identifier: BSD-3-Clause
// SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
/****************************************************************************
 *  dycelfem_config.cpp   -- see dycelfem_config.h
 ****************************************************************************/

#include "dycelfem_config.h"

#include <cstdlib>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>

using namespace std;

namespace
{
	map<string, string> g_values;
	vector<string>      g_order;      /* keys in file order, for dumpEffective */
	bool                g_loaded = false;
	string              g_path;
	string              g_dir;        /* directory containing the config */

	string trim(const string &s)
	{
		string::size_type b = s.find_first_not_of(" \t\r\n");
		if (b == string::npos) return "";
		string::size_type e = s.find_last_not_of(" \t\r\n");
		return s.substr(b, e - b + 1);
	}

	string lower(const string &s)
	{
		string r(s);
		for (string::size_type i = 0; i < r.size(); i++)
			r[i] = (char)tolower((unsigned char)r[i]);
		return r;
	}

	/* "view.autofit" -> "DYCELFEM_VIEW_AUTOFIT" */
	string envName(const string &key)
	{
		string r = "DYCELFEM_";
		for (string::size_type i = 0; i < key.size(); i++)
		{
			char c = key[i];
			r += (c == '.' || c == '-') ? '_' : (char)toupper((unsigned char)c);
		}
		return r;
	}

	const char *envLookup(const string &key)
	{
		const char *v = getenv(envName(key).c_str());
		if (v != NULL && v[0] != '\0') return v;
		return NULL;
	}

	string dirOf(const string &p)
	{
		string::size_type s = p.find_last_of("/\\");
		return (s == string::npos) ? string(".") : p.substr(0, s);
	}

	/* Shared by load() and looksLikeConfig(). */
	bool parseInto(const string &path, map<string, string> &out, vector<string> &order)
	{
		ifstream in(path.c_str());
		if (!in.good()) return false;

		string line, section;
		while (getline(in, line))
		{
			line = trim(line);
			if (line.empty() || line[0] == '#' || line[0] == ';') continue;

			if (line[0] == '[')
			{
				string::size_type e = line.find(']');
				if (e == string::npos) continue;
				section = lower(trim(line.substr(1, e - 1)));
				continue;
			}

			string::size_type eq = line.find('=');
			if (eq == string::npos) continue;

			string key = lower(trim(line.substr(0, eq)));
			string val = trim(line.substr(eq + 1));
			if (key.empty()) continue;

			/* strip a trailing inline comment, but not inside quotes */
			if (!val.empty() && val[0] != '"')
			{
				string::size_type h = val.find_first_of("#;");
				if (h != string::npos) val = trim(val.substr(0, h));
			}
			else if (val.size() >= 2 && val[0] == '"')
			{
				string::size_type q = val.find('"', 1);
				if (q != string::npos) val = val.substr(1, q - 1);
			}

			if (!section.empty()) key = section + "." + key;
			if (out.find(key) == out.end()) order.push_back(key);
			out[key] = val;
		}
		return true;
	}
}

namespace dycfg
{

bool load(const string &p)
{
	g_values.clear();
	g_order.clear();
	g_loaded = parseInto(p, g_values, g_order);
	if (g_loaded)
	{
		g_path = p;
		g_dir  = dirOf(p);
	}
	return g_loaded;
}

bool        loaded() { return g_loaded; }
std::string path()   { return g_path; }

bool looksLikeConfig(const string &p)
{
	map<string, string> m;
	vector<string>      o;
	if (!parseInto(p, m, o)) return false;
	return m.find("input") != m.end();
}

bool has(const string &key)
{
	if (envLookup(key) != NULL) return true;
	return g_values.find(key) != g_values.end();
}

string getStr(const string &key, const string &def)
{
	const char *e = envLookup(key);
	if (e != NULL) return string(e);
	map<string, string>::const_iterator i = g_values.find(key);
	return (i == g_values.end()) ? def : i->second;
}

int getInt(const string &key, int def)
{
	string v = getStr(key, "");
	if (v.empty()) return def;
	istringstream ss(v);
	int r;
	return (ss >> r) ? r : def;
}

double getNum(const string &key, double def)
{
	string v = getStr(key, "");
	if (v.empty()) return def;
	istringstream ss(v);
	double r;
	return (ss >> r) ? r : def;
}

bool getBool(const string &key, bool def)
{
	string v = lower(getStr(key, ""));
	if (v.empty()) return def;
	if (v == "1" || v == "true"  || v == "yes" || v == "on")  return true;
	if (v == "0" || v == "false" || v == "no"  || v == "off") return false;
	return def;
}

string resolve(const string &filename)
{
	if (filename.empty()) return filename;
	if (filename[0] == '/' || filename[0] == '\\') return filename;      /* absolute */
	if (filename.size() > 1 && filename[1] == ':') return filename;      /* C:\... */
	if (!g_loaded || g_dir.empty() || g_dir == ".") return filename;
	return g_dir + "/" + filename;
}

void dumpEffective()
{
	cout << "----------------------------------------------------------------" << endl;
	if (g_loaded) cout << "DyCelFEM configuration: " << g_path << endl;
	else          cout << "DyCelFEM configuration: (none - built-in defaults)" << endl;

	for (size_t i = 0; i < g_order.size(); i++)
	{
		const string &k = g_order[i];
		const char *e = envLookup(k);
		cout << "  " << k << " = " << (e ? string(e) : g_values[k]);
		if (e) cout << "   [overridden by " << envName(k) << "]";
		cout << endl;
	}
	cout << "----------------------------------------------------------------" << endl;
}

}

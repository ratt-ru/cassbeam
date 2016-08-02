/* cassbeam - a Cassegrain antenna simulator
    Copyright (C) August 18, 2003  Walter Brisken

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, see
    <http://www.gnu.org/licenses/>.
 */
#ifndef __KEYVALUE__
#define __KEYVALUE__

#define KV_MAXPARMS	256	/* Max number of parameters for a KeyValue */
#define KV_INTERR	-999999
#define KV_FLOATERR	-999999.0

#include "vector.h"

struct KeyValue
{
	char *key[KV_MAXPARMS];
	char *value[KV_MAXPARMS];
	int nparms;
};

struct KeyValue *newKeyValue();
void deleteKeyValue(struct KeyValue *p);
void KeyValueaddparm(struct KeyValue *p, const char *key, const char *value);
void KeyValueupdateparm(struct KeyValue *p, const char *key, const char *value);
void KeyValueupdateparmdouble(struct KeyValue *p, const char *key, double value);
struct KeyValue *loadKeyValue(const char *filename);
int saveKeyValue(const struct KeyValue *kv, const char *filename);
int KeyValuekeyindex(const struct KeyValue *p, const char *key);
int getKeyValueint(const struct KeyValue *p, const char *key);
double getKeyValuedouble(const struct KeyValue *p, const char *key);
Vector getKeyValueVector(const struct KeyValue *p, const char *key);
const char *getKeyValuestring(const struct KeyValue *p, const char *key);
void printKeyValue(const struct KeyValue *p);

#endif

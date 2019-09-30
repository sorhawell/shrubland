#ifndef _FOO_H 
#define _FOO_H
#endif

#include <Arduino.h>

void sprint(int x) {Serial.print(x);}
//void sprint(uint16_t x) {Serial.print(x);}
//void sprint(float x) {Serial.print(x);}
void sprint(const char x[]) {Serial.print(x);}
//void sprint(double x) {Serial.print(x);}


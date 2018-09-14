
#include<Servo.h>
Servo s[7];

void setup() {
  for(int i=0;i<7;++i)
    s[i].attach(2+i,300,2500);
  Serial.begin(9600);
}

void loop() {
  if(Serial.available())
  {
    int d=Serial.read()-'a';
    int ang=Serial.parseInt(),ang2=0;
    if(d==2)
    {
      ang2=2500-(ang-500);
      s[6].writeMicroseconds(ang2);
    }
    s[d].writeMicroseconds(ang);
    Serial.println(d);
    Serial.println(ang);
  }
}

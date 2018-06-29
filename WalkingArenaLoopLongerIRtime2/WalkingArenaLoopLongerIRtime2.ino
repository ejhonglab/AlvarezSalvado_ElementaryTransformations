/*
  WalkingArenaLoop
  Loop to run walking arena
  
  read pin 0 for light and valve control (high = on)
  read pin 1 for wind control (high = on)
  write pins 20-21 for light and valve control (LIGHT -> low = on)(SSR+VALVES -> low = off) 
  write pin 19 for wind control(high = on)
  strobe pin 11 for IR light (low = on)
  trigger pin 12 for camera (low = on)

 */
 
int input = 0;  // Pin 0 light in (and odour)
int windin = 1; //pin 1 wind in

int IRout = 11;  // Pin 11 IR out
int cameraout = 12;  // Pin 12 camera out


int lightout = 21;  // Pin 21  light out
int valveout = 20;  // Pin 20  valves relay out
int compout = 4; //separate control for the compensate valve
int windout = 19;   // Pin 19 wind out

// the setup routine runs once when you press reset:
void setup() {                
  // initialize the digital pin as an output.
  pinMode(input, INPUT);  
  pinMode(windin, INPUT);   
 
  pinMode(IRout, OUTPUT); 
  pinMode(cameraout, OUTPUT); 
  
  pinMode(lightout, OUTPUT); 
  pinMode(valveout, OUTPUT); 
  pinMode(compout, OUTPUT);
  pinMode(windout, OUTPUT);
  
  digitalWrite(IRout, HIGH); 
  digitalWrite(cameraout, HIGH); 
  digitalWrite(lightout, HIGH); 
  digitalWrite(valveout, LOW); 
  digitalWrite(compout, LOW);
  digitalWrite(windout, LOW);
}

// the loop routine runs over and over again forever:
void loop() {
  // read and write values for light and valves
  int lightval = digitalRead(input);
  if (lightval==HIGH)
    {digitalWrite(lightout,LOW);}
  else if(lightval==LOW)
    {digitalWrite(lightout,HIGH);}
    
int valveval = digitalRead(input);
int windval = digitalRead(windin);
Serial.print(windval);
if (windval==HIGH&&valveval==HIGH)
  {digitalWrite(windout,HIGH);
  digitalWrite(compout,LOW);
  digitalWrite(valveout,HIGH);}
  else if (windval==HIGH)
  {digitalWrite(windout,HIGH);
  digitalWrite(compout,HIGH);
  digitalWrite(valveout,LOW);}
  else
  {digitalWrite(windout,LOW);
  digitalWrite(compout,LOW);
  digitalWrite(valveout,LOW);}
  

//    
//int valve2val = digitalRead(valve2in);
//  if (valve2val==HIGH)
//    {digitalWrite(valve2out,HIGH);}
//  else if(valve2val==LOW)
//    {digitalWrite(valve2out,LOW);}
//    
//int valve3val = digitalRead(valve3in);
//  if (valve3val==HIGH)
//    {digitalWrite(valve3out,HIGH);}
//  else if(valve3val==LOW)
//    {digitalWrite(valve3out,LOW);}
    
  // trigger IR light and camera  
  digitalWrite(IRout,LOW);
  digitalWrite(cameraout,LOW);
  delay(1);
  
  // resets camera
  digitalWrite(cameraout,HIGH);
  delay(14);

  // reset IR light
  digitalWrite(IRout,HIGH);

  // wait 10ms (runs at 50Hz)
  delay(5);
}








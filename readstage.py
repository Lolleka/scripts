#!/usr/bin/python
import serial as s
import curses

ser = s.Serial(port='/dev/ttyUSB0',baudrate=115200)

ser.write("0 1 setcdfunc\r\n".encode())
ser.write("0 2 setcdfunc\r\n".encode())

ser.write("2 np\r\n".encode())	
l = ser.readline().decode()[:-2]
print(l)
ser.write("1 np\r\n".encode())
p = ser.readline().decode()[:-2]
print(p)

ser.close()
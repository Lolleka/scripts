import serial

ser = serial.Serial("/dev/ttyS10",9600)

ser.close()
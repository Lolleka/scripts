#!/usr/bin/python
import serial as s
import curses


class Position:
    def __init__(self,x,y):
        self.x = x
        self.y = y

ser = s.Serial(port='/dev/ttyUSB0',baudrate=115200)

ser.write("0 1 setcdfunc\r\n".encode())
ser.write("0 2 setcdfunc\r\n".encode())

stdscr = curses.initscr()

height,width = stdscr.getmaxyx()

curses.cbreak()
curses.noecho()

stdscr.keypad(1)

stdscr.addstr(0, 0, "Hit 'q' to quit, 's' to enter save mode, 'l' to enter load mode")
stdscr.addstr(3, 20, "Up")
stdscr.addstr(5, 20, "Down")
stdscr.addstr(4, 10, "Left")
stdscr.addstr(4, 30, "Right")


key=''
inc = 1.0
stdscr.addstr(2, 0, "Increment: " + str(inc))
stdscr.refresh()

cpos = Position(0,0) #Current Position

spos = [None]*8

for i in range(0,8):
	spos[i] = Position(0,0) #Saved Position


while key != ord('q'):
	key = stdscr.getch()
	#stdscr.addch(20,25,key)
#	stdscr.refresh()
	if key == curses.KEY_UP:
		stdscr.addstr(3, 20, "Up")
		ser.write((str(inc) + " 2 nr\r\n").encode())
	elif key == curses.KEY_DOWN:
		stdscr.addstr(5, 20, "Down")
		ser.write(("-" + str(inc) + " 2 nr\r\n").encode())		
	elif key == curses.KEY_LEFT:
		stdscr.addstr(4, 10, "Left")
		ser.write(("-" + str(inc) + " 1 nr\r\n").encode())
	elif key == curses.KEY_RIGHT:
		stdscr.addstr(4, 30, "Right")
		ser.write((str(inc) + " 1 nr\r\n").encode())
	elif key == ord('a'):
		inc = inc * 2
		stdscr.addstr(2, 1, "inc: " + str(inc))
	elif key == ord('z'):
		inc = inc / 2
		stdscr.addstr(2, 1, "inc: " + str(inc))
	elif key == ord('s'):
		stdscr.addstr(25,0,"Choose save slot [1-8]: ")
		stdscr.refresh()	
		key2 = stdscr.getch()
		if key2 >= ord('1') and key2 <= ord('8'):
			#Read stage position
			ser.write("1 np\r\n".encode())
			cpos.x = float(ser.readline().decode()[:-2])
			ser.write("2 np\r\n".encode())
			cpos.y = float(ser.readline().decode()[:-2])
			index = key2 - ord('1')
			spos[index].x = cpos.x
			spos[index].y = cpos.y
		stdscr.move(25,0)
		stdscr.clrtoeol()
		stdscr.refresh()
	elif key == ord('l'):
		stdscr.addstr(height-1,0,"Choose slot [1-8]: ")
		stdscr.refresh()
		key2 = stdscr.getch()
		if key2 >= ord('1') and key2 <= ord('8'):
			index = key2 - ord('1')
			ser.write((str(spos[index].x) + " 1 nm\r\n").encode())
			ser.write((str(spos[index].y) + " 2 nm\r\n").encode())
		stdscr.move(height-1,0)
		stdscr.clrtoeol()
		stdscr.refresh()

	#Read stage position
	ser.write("1 np\r\n".encode())
	cpos.x = float(ser.readline().decode()[:-2])
	ser.write("2 np\r\n".encode())
	cpos.y = float(ser.readline().decode()[:-2])
	stdscr.addstr(8, 1, "Current Position: X = " + str(cpos.x) + " : Y = " + str(cpos.y))
	for i in range(0,8):
		stdscr.addstr(9+i,1, "Saved " + str(i+1) + ": X = " + str(spos[i].x) + " : Y = " + str(spos[i].y))

	stdscr.refresh()

curses.endwin()
ser.close()
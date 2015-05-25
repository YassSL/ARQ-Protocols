#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream>
#include <math.h> 
#include <vector>
#include <utility>      // std::pair, std::make_pair
#include <queue>          // std::priority_queue
#include <fstream>



using namespace std;

double bit_generator(double p) {
    double y;
 	y = ((double)rand()/(double)RAND_MAX);			
	 
	if ( y >= p) { 
		return 1;
	} else {
		return 0;
	}
}


enum event_type {Ack = 0, time_out};


class Event { 
	event_type type; 
	double time; 
	bool errorFlag;
	bool lossFlag;
	int SN; 
public:
	Event(event_type, double, bool, bool, int);
	Event(event_type, double);
	event_type getEventType() { return type; }
	void setEventType(event_type eType) { type = eType; } 
	double getTime() { return time; }
	void setTime(int eTime) { time = eTime; }
	bool getErrorFlag() { return errorFlag; } 
	void setErrorFlag(bool flag) { errorFlag = flag; } 
	bool getLossFlag() { return lossFlag; } 
	void setLossFlag(bool flag) { lossFlag = flag; }
	int getSNumber() {return SN; }
	void setSNumber(int number) {SN = number; }
};

Event::Event(event_type eType, double eTime, bool eFlag, bool lFlag, int sNumber) {
	type = eType;
	time = eTime;
	errorFlag = eFlag;
	lossFlag = lFlag;
	SN = sNumber;
};

Event::Event(event_type eType, double eTime) { 
	type = eType;
	time = eTime;
};

struct sort_func {
    bool operator()(const pair<Event, double> &left, const pair<Event, double> &right) {
        return left.second > right.second;
    }
};




priority_queue<pair<Event, double>, vector<pair<Event, double> >, sort_func > DES;

static const int deliveredPackets = 10000;
static const double l = 12000.0;	// Average packet length in bits
static const double H = 432.0;	// Average packet length in bits
static const double L = l + H;	// Average packet length in bits
static const double C = 5000000.0;	// Transmission rate of the output link in bits/sec

static double tp = 0.005;	// Channel Propagation Delay
static double delta_t = 2.5 * tp;
static double BER = 0; // bit error rate
static double tc = 0;	// Current time

static int SN = 0;
static int Next_Expected_Ack = 1;
static int Next_Expected_Frame = 0;
static int numOfDeliveredPackets = 0;

Event channel(double currTime, int SN, double packetLen) {
	int zeroCounter = 0;
	int bit;
	Event ch(Ack, currTime + packetLen / C + tp, false, false, SN);
	for (int i = 0; i < packetLen; i++) {
		bit = bit_generator(BER);
		if (bit == 0) {
			zeroCounter++; 
		}
	} 
		
	if(zeroCounter >= 5) {
		ch.setLossFlag(true);
	} else if (zeroCounter < 5 && zeroCounter > 0) { 
		ch.setErrorFlag(true);
	}
	
	return ch;
}

int reciever(Event chOut) { 
	if(!chOut.getErrorFlag()) {
		if(chOut.getSNumber() == Next_Expected_Frame) { 
			Next_Expected_Frame = (Next_Expected_Frame + 1) % 2;
		}
	}
	return Next_Expected_Frame;

}

Event send(double currTime, int SN, double packet_len) {
	Event forwardChannel = channel(currTime, SN, packet_len);
	int RN = reciever(forwardChannel);
	Event reverseChannel = channel(forwardChannel.getTime(), RN, H);
	bool finalErrorFlag = forwardChannel.getErrorFlag() || reverseChannel.getErrorFlag(); 
	bool finalLossFlag = forwardChannel.getLossFlag() || reverseChannel.getLossFlag();
	Event finalAck(Ack, reverseChannel.getTime(), finalErrorFlag, finalLossFlag, RN);
	return finalAck;
} 

double sender() { 
	while(numOfDeliveredPackets < deliveredPackets) {
		while(!DES.empty()) {
			DES.pop();
		}
		Event t_out(time_out, tc + L/C + delta_t);
		DES.push(make_pair<Event, double>(t_out, t_out.getTime()));	 // Register time_out event
	
		Event finalAck = send(tc, SN, L);
	
		if (!finalAck.getLossFlag() {
			DES.push(make_pair<Event, double>(finalAck, finalAck.getTime()));
		}
	
		pair<Event, double> event_pair = DES.top(); // Get the first element
		Event currEvent = event_pair.first;
		tc = event_pair.second;
		
		DES.pop();
		if(currEvent.getEventType() != time_out) {	// ACk was dequeued
			if (!currEvent.getErrorFlag() && currEvent.getSNumber() == Next_Expected_Ack) {
				numOfDeliveredPackets++;
				Next_Expected_Ack = (Next_Expected_Ack + 1) % 2;
				SN = (SN + 1) % 2;
			}  else { // if ACK recieved was in error, get the next event which is time_out and update the current time
				event_pair = DES.top();
				if (event_pair.first.getEventType() == time_out) {
					tc = event_pair.second;
				}
				DES.pop();
			}
		}
	
	}
	return numOfDeliveredPackets * l / tc;
}



int main(int argc, char* argv[]) { 
	srand(time(0));
	double bers[3] = {0, 0.00001, 0.0001};
	double tps[2] = {0.005, 0.25};
	ofstream ABPFILE("ABP.csv"); 
	for (double i = 2.5; i < 15; i= i +2.5) {
		for (int k = 0; k < 2; k++) {
			for(int j = 0; j < 3; j++) {
				BER = bers[j];
				tp = tps[k];
				delta_t = i * tp;
				tc = 0;	// Current time
				SN = 0;
			    Next_Expected_Ack = 1;
				Next_Expected_Frame = 0;
				numOfDeliveredPackets = 0;
				ABPFILE << sender();
				ABPFILE<< ",";

			}
		}
		ABPFILE << "\n";
	}
	ABPFILE.close();
	return 0;	
}

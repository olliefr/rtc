-- Send rotor rig firmware .bin file via xmodem and run it
-- Oliver Frolovs, University of Bristol, 2019

lua_senddata("loady 0x80000000; go 0x80000000", true);
lua_ymodem_snd("C:\\Local\\rtc\\rigs\\rotor\\rotor.bin");

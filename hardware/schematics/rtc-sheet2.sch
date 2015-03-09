v 20130925 2
C 40000 40000 0 0 0 title-B.sym
C 40900 40700 1 0 0 header46.sym
{
T 41200 50000 5 10 1 1 0 0 1
refdes=J2
T 41600 50000 5 10 1 1 0 0 1
device=P8
T 40900 40700 5 10 0 1 0 0 1
footprint=HEADER46_2
}
C 43800 40700 1 0 0 header46.sym
{
T 44100 50000 5 10 1 1 0 0 1
refdes=J3
T 44500 50000 5 10 1 1 0 0 1
device=P9
T 43800 40700 5 10 0 1 0 0 1
footprint=HEADER46_2
}
C 48500 47500 1 0 0 ad5064-1.sym
{
T 50400 49700 5 10 1 1 0 6 1
refdes=U3
T 49500 48800 5 10 0 0 0 0 1
device=AD5064-1
T 49500 49000 5 10 0 0 0 0 1
footprint=TSSOP14
}
C 54300 48800 1 0 1 adr4540.sym
{
T 54000 49900 5 10 1 1 0 6 1
refdes=U4
T 53900 49600 5 10 1 1 0 6 1
device=ADR4540
T 51700 48600 5 10 0 1 0 6 1
footprint=SO8
}
C 52700 48500 1 90 0 capacitor-1.sym
{
T 52000 48700 5 10 0 0 90 0 1
device=CAPACITOR
T 52300 48700 5 10 1 1 90 0 1
refdes=C11
T 51800 48700 5 10 0 0 90 0 1
symversion=0.1
T 52700 48500 5 10 0 0 90 0 1
value=100nF
T 52700 48500 5 10 0 1 0 0 1
footprint=0805
}
C 53400 48200 1 0 0 gnd-1.sym
C 49500 50300 1 0 0 vcc-1.sym
C 54300 49400 1 0 0 vcc-1.sym
C 54700 48500 1 90 0 capacitor-1.sym
{
T 54000 48700 5 10 0 0 90 0 1
device=CAPACITOR
T 54300 48700 5 10 1 1 90 0 1
refdes=C12
T 53800 48700 5 10 0 0 90 0 1
symversion=0.1
T 54700 48500 5 10 0 0 90 0 1
value=100nF
T 54700 48500 5 10 0 1 0 0 1
footprint=0805
}
N 54300 49400 54500 49400 4
N 52500 48500 54500 48500 4
N 53500 48800 53500 48500 4
N 50800 49400 52700 49400 4
C 49700 49900 1 0 0 capacitor-1.sym
{
T 49900 50600 5 10 0 0 0 0 1
device=CAPACITOR
T 49900 50300 5 10 1 1 0 0 1
refdes=C13
T 49900 50800 5 10 0 0 0 0 1
symversion=0.1
T 49700 49900 5 10 0 0 0 0 1
value=100nF
T 49700 49900 5 10 0 0 0 0 1
footprint=0805
}
C 50700 49800 1 0 0 gnd-1.sym
C 49600 47100 1 0 0 gnd-1.sym
N 49700 47400 49700 47600 4
N 50800 50100 50600 50100 4
N 49700 50300 49700 49900 4
C 47700 47800 1 90 0 vcc-1.sym
N 47700 48000 48600 48000 4
C 45800 48300 1 270 1 vcc-1.sym
C 43600 48300 1 90 0 vcc-1.sym
C 43600 49100 1 90 0 vdd-1.sym
C 45400 49100 1 270 1 vdd-1.sym
C 45700 49600 1 90 0 gnd-1.sym
C 43300 49600 1 270 1 gnd-1.sym
C 45700 41000 1 90 0 gnd-1.sym
C 43300 41000 1 270 1 gnd-1.sym
N 43800 40900 43800 41300 4
N 43800 41100 43600 41100 4
N 45200 41300 45200 40900 4
N 45200 41100 45400 41100 4
N 43600 48500 43800 48500 4
N 45800 48500 45200 48500 4
N 43600 49300 43800 49300 4
N 45200 49300 45400 49300 4
N 45400 49700 45200 49700 4
N 43600 49700 43800 49700 4
C 54500 45300 1 0 0 header12.sym
{
T 54800 47800 5 10 1 1 0 0 1
refdes=J4
T 54500 45300 5 10 0 0 0 0 1
footprint=HEADER12_2
}
C 54100 47700 1 0 0 vcc-1.sym
N 54500 47100 52000 47100 4
N 52000 47100 52000 49400 4
N 54500 47500 54300 47500 4
N 54300 47500 54300 47700 4
C 56500 46400 1 90 0 gnd-1.sym
N 55900 47500 55900 45500 4
N 56200 46500 55900 46500 4
N 40300 46500 40900 46500 4
{
T 40300 46500 5 10 1 1 0 0 1
netname=OS1
}
N 40300 47700 40900 47700 4
{
T 40300 47700 5 10 1 1 0 0 1
netname=BUSY
}
N 40900 46900 40300 46900 4
{
T 40300 46900 5 10 1 1 0 0 1
netname=STBY
}
N 42300 46500 42900 46500 4
{
T 42500 46500 5 10 1 1 0 0 1
netname=OS0
}
N 42900 46900 42300 46900 4
{
T 42500 46900 5 10 1 1 0 0 1
netname=OS2
}
N 42300 47300 42900 47300 4
{
T 42300 47300 5 10 1 1 0 0 1
netname=RANGE
}
N 42900 47700 42300 47700 4
{
T 42300 47700 5 10 1 1 0 0 1
netname=RESET
}
N 42300 48100 42900 48100 4
{
T 42200 48100 5 10 1 1 0 0 1
netname=CONVST
}
N 43200 43700 43800 43700 4
{
T 43200 43700 5 10 1 1 0 0 1
netname=SCLK1
}
N 45800 44500 45200 44500 4
{
T 45500 44500 5 10 1 1 0 0 1
netname=CS1
}
N 45800 44100 45200 44100 4
{
T 45300 44100 5 10 1 1 0 0 1
netname=MISO1
}
N 42900 48500 42300 48500 4
{
T 42600 48500 5 10 1 1 0 0 1
netname=CLR
}
N 40300 48500 40900 48500 4
{
T 40300 48500 5 10 1 1 0 0 1
netname=LDAC
}
N 45800 45700 45200 45700 4
{
T 45200 45700 5 10 1 1 0 0 1
netname=SCLK0
}
N 43200 46500 43800 46500 4
{
T 43200 46500 5 10 1 1 0 0 1
netname=CS0
}
N 43200 45700 43800 45700 4
{
T 43200 45700 5 10 1 1 0 0 1
netname=MOSI0
}
C 45400 48500 1 270 0 capacitor-2.sym
{
T 46100 48300 5 10 0 0 270 0 1
device=POLARIZED_CAPACITOR
T 45800 48300 5 10 1 1 270 0 1
refdes=C4
T 46300 48300 5 10 0 0 270 0 1
symversion=0.1
T 45400 48500 5 10 1 1 270 0 1
value=10uF
T 45400 48500 5 10 0 1 0 0 1
footprint=0805
}
C 45500 47100 1 0 0 gnd-1.sym
N 45600 47400 45600 47600 4
N 50800 49100 51400 49100 4
{
T 50800 49100 5 10 1 1 0 0 1
netname=VOUTA
}
N 50800 48900 51400 48900 4
{
T 50800 48900 5 10 1 1 0 0 1
netname=VOUTB
}
N 50800 48700 51400 48700 4
{
T 50800 48700 5 10 1 1 0 0 1
netname=VOUTC
}
N 50800 48500 51400 48500 4
{
T 50800 48500 5 10 1 1 0 0 1
netname=VOUTD
}
N 54500 45900 53800 45900 4
{
T 53800 45900 5 10 1 1 0 0 1
netname=VOUTD
}
N 54500 45500 53800 45500 4
{
T 53800 45500 5 10 1 1 0 0 1
netname=VOUTB
}
N 48000 49100 48600 49100 4
{
T 48000 49100 5 10 1 1 0 0 1
netname=CLR
}
N 48000 49300 48600 49300 4
{
T 48000 49300 5 10 1 1 0 0 1
netname=LDAC
}
N 48000 48400 48600 48400 4
{
T 48000 48400 5 10 1 1 0 0 1
netname=CS0
}
N 48000 48200 48600 48200 4
{
T 48000 48200 5 10 1 1 0 0 1
netname=SCLK0
}
N 48000 48600 48600 48600 4
{
T 48000 48600 5 10 1 1 0 0 1
netname=MOSI0
}
C 52700 42200 1 0 0 ad8639.sym
{
T 52695 42200 5 10 0 1 0 0 1
footprint=SO8
T 52995 43600 5 10 1 1 0 0 1
device=AD8639
T 53800 43600 5 10 1 1 0 0 1
refdes=U6
}
C 46600 42200 1 0 0 icl7660-1.sym
{
T 48400 43700 5 10 0 0 0 0 1
device=icl7660
T 46600 42200 5 10 0 0 0 0 1
footprint=SO8
T 48300 43800 5 10 1 1 0 6 1
refdes=U5
}
C 49000 41500 1 90 0 capacitor-2.sym
{
T 48300 41700 5 10 0 0 90 0 1
device=POLARIZED_CAPACITOR
T 48100 41700 5 10 0 0 90 0 1
symversion=0.1
T 49000 41500 5 10 0 1 0 0 1
footprint=0805
T 48500 41700 5 10 1 1 90 0 1
refdes=C14
T 49200 41900 5 10 1 1 90 0 1
value=10uF
}
C 49300 43400 1 270 0 capacitor-2.sym
{
T 50000 43200 5 10 0 0 270 0 1
device=POLARIZED_CAPACITOR
T 50200 43200 5 10 0 0 270 0 1
symversion=0.1
T 49300 43400 5 10 0 1 0 0 1
footprint=0805
T 49800 43400 5 10 1 1 270 0 1
refdes=C16
T 49800 42900 5 10 1 1 270 0 1
value=10uF
}
C 49000 43500 1 90 0 capacitor-2.sym
{
T 48300 43700 5 10 0 0 90 0 1
device=POLARIZED_CAPACITOR
T 48100 43700 5 10 0 0 90 0 1
symversion=0.1
T 49000 43500 5 10 0 1 0 0 1
footprint=0805
T 48500 43700 5 10 1 1 90 0 1
refdes=C15
T 49100 44100 5 10 1 1 90 0 1
value=100nF
}
C 50800 42900 1 0 0 resistor-1.sym
{
T 51100 43300 5 10 0 0 0 0 1
device=RESISTOR
T 50800 42900 5 10 0 1 0 0 1
footprint=0603
T 51000 43200 5 10 1 1 0 0 1
refdes=R11
T 51300 42700 5 10 1 1 0 0 1
value=10k
}
N 51700 43000 52700 43000 4
N 51700 43000 51700 43300 4
N 55500 43000 55500 42700 4
N 55500 42700 54300 42700 4
N 54600 43000 54300 43000 4
N 52600 43300 52700 43300 4
N 49100 43400 49500 43400 4
N 49500 42500 49100 42500 4
N 49100 42500 49100 42800 4
N 49100 42800 48600 42800 4
N 48600 43100 49100 43100 4
N 49100 43100 49100 43400 4
C 46300 42200 1 0 0 gnd-1.sym
C 48700 41100 1 0 0 gnd-1.sym
C 48900 44800 1 180 0 gnd-1.sym
N 46400 42500 46600 42500 4
N 48800 41500 48800 41400 4
N 48600 42500 48800 42500 4
N 48800 42500 48800 42400 4
N 48800 44500 48800 44400 4
N 48600 43400 48800 43400 4
N 48800 43400 48800 43500 4
N 48800 42400 52700 42400 4
C 49100 43700 1 0 0 vcc-1.sym
C 54200 43300 1 0 0 vcc-1.sym
N 49300 43700 49300 43500 4
N 49300 43500 48800 43500 4
N 54300 43300 54400 43300 4
C 50600 43000 1 0 0 vss-1.sym
C 56200 42700 1 0 0 vss-1.sym
C 51700 43200 1 0 0 resistor-1.sym
{
T 52000 43600 5 10 0 0 0 0 1
device=RESISTOR
T 51700 43200 5 10 0 1 0 0 1
footprint=0603
T 51800 43500 5 10 1 1 0 0 1
refdes=R12
T 52300 43500 5 10 1 1 0 0 1
value=10k
}
C 54600 42900 1 0 0 resistor-1.sym
{
T 54900 43300 5 10 0 0 0 0 1
device=RESISTOR
T 54600 42900 5 10 0 1 0 0 1
footprint=0603
T 54700 43200 5 10 1 1 0 0 1
refdes=R13
T 55200 43200 5 10 1 1 0 0 1
value=10k
}
C 55500 42600 1 0 0 resistor-1.sym
{
T 55800 43000 5 10 0 0 0 0 1
device=RESISTOR
T 55500 42600 5 10 0 1 0 0 1
footprint=0603
T 55700 42900 5 10 1 1 0 0 1
refdes=R14
T 56000 42400 5 10 1 1 0 0 1
value=10k
}
N 54500 46700 52700 46700 4
N 52700 43300 52700 46700 4
N 54500 46300 53100 46300 4
N 53100 44200 53100 46300 4
N 53100 44200 54600 44200 4
N 54600 44200 54600 43000 4
N 52700 42700 52000 42700 4
{
T 52000 42700 5 10 1 1 0 0 1
netname=VOUTA
}
N 54300 42400 55000 42400 4
{
T 54300 42400 5 10 1 1 0 0 1
netname=VOUTC
}

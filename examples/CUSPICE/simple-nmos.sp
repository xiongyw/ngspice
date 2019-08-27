.title simple MOSFET circuit

* circuit is based on <http://web.engr.oregonstate.edu/~traylor/ece112/beamer_lectures/mosfet_spice.pdf>,
* with change of nmos model and modifications of vdd, v_input and rdrain.

.include "Circuits/technology/45nm_HP.pm"   ; BSIM4v7 model `nmos`

.options badchr=1 ingold=1 numdgt=4

v_input vin    gnd    0.0  pulse(0 1 10m 10m 10m 10m 50m)
Vdd     vdd    gnd    1   ; nominal Vdd of the nmos is 1V
Mq1     drain  vin    gnd  gnd   nmos  ; d g s b model
rdrain  vdd    drain  40k ; resistor in drain lead
.control
tran 1.0u 50m
plot v(vin) v(drain) xlimit 0 50m
.endc
.end

Netflix
=======

About
-----

This is a project I wrote as part of my undergraduate degree in Engineering at Cambridge (supervised by Zoubin Ghahramani).

The project was around using a Bayesian approach to enter the Netflix Competition (http://www.netflixprize.com/), 
a large collaborative filtering challenge in the field of machine learning. I was investigating different approaches,
including different bits of the dataset or not, and the project culminated in a dissertation. This is the
code I wrote to obtain my results.

In short, Netflix give a large number of ratings given by users for a set of movies. The challenge was to predict other
genuine ratings on Netflix using that data, and to improve on the predictions of Netflix's own algorithm by 10%. The best
submission I got in was roughly a 5.2% improvement on Netflix's own algorithm (not bad!)

Cool Feature
------------

One cool feature is the innovative way I packed bits in. A lot of metadata about each reading is bundled into a single
64-bit long. This is done by having a set of distinct prime "keys" which are used as the right-hand side of the modulo
operator to get the data out of the long. So if the key for favourite_number is 11 and the key for birthday_month is 13
and the 64-bit long is 97, the favourite_number value is 9 (97 % 11) and the birthday_month value is 6 (97 % 13). The 
only tricky bit is creating the 64-bit long from the values it has to hold: that's done by NetflixDataNG::GetBigRating.

Data
----

I did some pre-processing on the dataset (which I think you can still download) before reading it into this program. This is
contained within PreProc.cpp and .h.

Compilation
-----------

	cd Netflix
	mkdir build
	cd build
	cmake ..
	make

Credits
-------

Zoubin Ghahramani and Sinead Williamson both helped in supervising me. The code began as a C++ port of Vibes by John Winn 
(http://johnwinn.org/). It was ported because Java isn't efficient enough to fit the whole Netflix dataset in memory at once,
whereas C++ is on my 8GB RAM desktop. Where Vibes is a general code with a GUI, this code is designed for the single purpose
of using variational inference to solve a particular model for the Netflix problem.

Disclaimer
----------

I'm afraid this is largely PhD-ware. If you have any questions, I'm happy to answer them.

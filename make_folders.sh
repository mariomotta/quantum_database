rm -rf auxiliary_material first_quantization second_quantization

mkdir first_quantization
	mkdir first_quantization/pad
		mkdir first_quantization/pad/operators
		mkdir first_quantization/pad/projection_after_variation
			mkdir first_quantization/pad/projection_after_variation/circuits
				mkdir first_quantization/pad/projection_after_variation/circuits/cascade_3
				mkdir first_quantization/pad/projection_after_variation/circuits/cascade_4
				mkdir first_quantization/pad/projection_after_variation/circuits/cascade_5

	        mkdir first_quantization/pad/variation_after_projection
		        mkdir first_quantization/pad/variation_after_projection/circuits
	                        mkdir first_quantization/pad/variation_after_projection/circuits/cascade_3
                                mkdir first_quantization/pad/variation_after_projection/circuits/cascade_4
                                mkdir first_quantization/pad/variation_after_projection/circuits/cascade_5
				# aggiungere l'Ry
				# aggiungere il qUCCSD, automaticamente ristretto, ds

	mkdir first_quantization/trim
	        mkdir first_quantization/trim/operators
		mkdir first_quantization/trim/circuits
			mkdir first_quantization/trim/circuits/cascade_3
                        mkdir first_quantization/trim/circuits/cascade_4
                        mkdir first_quantization/trim/circuits/cascade_5

# ------------------------

mkdir second_quantization
	mkdir second_quantization/operators
	mkdir second_quantization/circuits
		mkdir second_quantization/circuits/hardware_efficient
			mkdir second_quantization/circuits/hardware_efficient/ry_linear_1
                        mkdir second_quantization/circuits/hardware_efficient/ry_linear_2
                        mkdir second_quantization/circuits/hardware_efficient/ry_linear_3
                        mkdir second_quantization/circuits/hardware_efficient/ry_linear_4
                        mkdir second_quantization/circuits/hardware_efficient/ry_linear_5
                        mkdir second_quantization/circuits/hardware_efficient/ry_linear_6
                        mkdir second_quantization/circuits/hardware_efficient/ry_linear_7
                        mkdir second_quantization/circuits/hardware_efficient/ry_linear_8
                        mkdir second_quantization/circuits/hardware_efficient/ry_linear_9

                        mkdir second_quantization/circuits/hardware_efficient/ry_full_1
                        mkdir second_quantization/circuits/hardware_efficient/ry_full_2
                        mkdir second_quantization/circuits/hardware_efficient/ry_full_3
                        mkdir second_quantization/circuits/hardware_efficient/ry_full_4
                        mkdir second_quantization/circuits/hardware_efficient/ry_full_5
                        mkdir second_quantization/circuits/hardware_efficient/ry_full_6
                        mkdir second_quantization/circuits/hardware_efficient/ry_full_7
                        mkdir second_quantization/circuits/hardware_efficient/ry_full_8
                        mkdir second_quantization/circuits/hardware_efficient/ry_full_9

		mkdir second_quantization/circuits/qUCCSD
	                mkdir second_quantization/circuits/qUCCSD/unrestricted
				mkdir second_quantization/circuits/qUCCSD/unrestricted/singles_doubles
                                mkdir second_quantization/circuits/qUCCSD/unrestricted/doubles_singles
                                # aggiungere l'ordinamento automatico
	                mkdir second_quantization/circuits/qUCCSD/restricted
                                mkdir second_quantization/circuits/qUCCSD/restricted/singles_doubles
                                mkdir second_quantization/circuits/qUCCSD/restricted/doubles_singles
				# aggiungere l'ordinamento automatico

# ------------------------

mkdir auxiliary_material
	mkdir auxiliary_material/scripts
	mkdir auxiliary_material/scf_fci # +
	mkdir auxiliary_material/figures
	mkdir auxiliary_material/writeup

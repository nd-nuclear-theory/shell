          // store matrix element
          //
          // H2 external storage is as NAS, but internal storage is as
          // AS.  So normalization conversion on input is "NASToAS".
          double conversion_factor = 1.;
          if (bra.index1()==bra.index2())
            conversion_factor *= (sqrt(2.));
          if (ket.index1()==ket.index2())
            conversion_factor *= (sqrt(2.));
          double matrix_element = conversion_factor * input_matrix_element;


          // retrieve matrix element for output
          //
          // H2 external storage is as NAS, but internal storage is as
          // AS.  So normalization conversion on output is "ASToNAS".
          double conversion_factor = 1.;
          if (bra.index1()==bra.index2())
            conversion_factor *= 1/(sqrt(2.));
          if (ket.index1()==ket.index2())
            conversion_factor *= 1/(sqrt(2.));
          double matrix_element = matrix(bra_index,ket_index);
          float output_matrix_element = conversion_factor * matrix_element;


package com.arbonaut.sampling;


/**
 * Exception class thrown if supplied geometry has missing or incompatible 
 * coordinate spatial reference. 
 */
public class BadCRSException extends Exception {
    
    public BadCRSException(String message) {
        super(message);
    }
}
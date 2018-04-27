
import turtle
from getStruct import*
from hivSeqs import*

def render(RNA, matches):
    """draws the RNA sequence graphically, showing every nucleotide in the RNA string
    and showing the matches in the positions given in the list matches."""
    turtle.dot(5, "blue")
    turtle.write(RNA[0], font = ("arial", 18, "normal"))
    for i in range(len(RNA)):
        if i != 0:
            turtle.forward(15)
            turtle.dot(5, "blue")
            turtle.write(RNA[i], font = ("arial", 18, "normal"))
    turtle.penup()
    turtle.setx(0)
    turtle.pendown()
    turtle.pencolor("red")
    for i in range(len(matches)):
        turtle.penup()
        turtle.setx(matches[i][0] * 15)
        turtle.pendown()
        turtle.right(90)
        turtle.forward(0.5 * 15 * (matches[i][1]-matches[i][0]))
        turtle.left(90)
        turtle.forward(15 * (matches[i][1]-matches[i][0]))
        turtle.left(90)
        turtle.forward(0.5 * 15 * (matches[i][1]-matches[i][0]))
        turtle.right(90)
    turtle.penup()
    turtle.setx(0)
    turtle.pendown()
        
def showStruct(RNA):
    """takes an RNA string as input, finds an optimal folding (by calling the getStruct
    function that you wrote in the previous problem), and finally rendering that fold."""
    memo = {}
    coolOption = getStruct(RNA, memo)
    render(RNA, coolOption[1])
    
#All four s1RegionSamples had broadly similar structure, adding additional support to the
#hypothesis that the pairings in this region have a significant biological role